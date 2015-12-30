function [ S_est, Img_est ] = general_solver( decorr,sampling_mecanism,method,SNR,Img,H,sources, N,n1,n2,I,J,nb_meas )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



%% Compressed sampling step

    if (decorr==1) && (~strcmp(sampling_mecanism ,'Block_diag') )
        error('Error: Decorelation applies only for Block_diag sampling')
    end  

    %Sampling operator identification.
    if strcmp(sampling_mecanism ,'Block_diag')
       opCS = opRC_measur_matrix( N, nb_meas);   % 2D Image CS by Random Convolution RC(J. Romberg 2009).
       opCS_multichannel = opBlock_diag_same(J,opCS );
          % t_gen = sigma_gen(N);
          % s_gen = theta(N);
          % opCS_multichannel = opBlock_Diag_RC_measur_matrix(N,J,nb_meas,s_gen, t_gen);
        Phi = @(x) opCS_multichannel(x,1);
        Phit = @(x) opCS_multichannel(x,2);
        nu_phi = 1;    % Since opRC is a normalized tight frame.
        tight_phi = 1; % Phi is a tight frame
        
        if decorr == 1 %    if decorrelation step applies we need to redefine
                       %    the block diagonal operator Phi with a smaller
                       %    dimension.

            %opCS_multichannel1 = opBlock_Diag_RC_measur_matrix(N,I,nb_meas,s_gen, t_gen);
            opCS_multichannel1 = opBlock_diag_same(I,opCS );
            Phi1 = @(x) opCS_multichannel1(x,1);
            Phit1 = @(x) opCS_multichannel1(x,2);
        else
            Phi1 = Phi;     % Here no difference between Phi1 and Phi. This
                            % is done for compatibilty reason
            Phit1 = Phit;
        end


    elseif strcmp(sampling_mecanism ,'DBD')

        opCS_multichannel = opDBD_RC_measur_matrix(N,J , nb_meas); % "Distributed independent" sampling by RC
        Phi = @(x) opCS_multichannel(x,1);
        Phit = @(x) opCS_multichannel(x,2);
        nu_phi = 1;         % Since opRC is a normalized tight frame.
        tight_phi = 1;      % Phi is a tight frame
        Phi1 = Phi;         % Here no difference between Phi1 and Phi. This
                            % is done for compatibilty reason
        Phit1 = Phit;

    elseif strcmp(sampling_mecanism ,'Dense')   
        opCS = opRC_measur_matrix( N*J, nb_meas*J);   % Dense CS by RC
        Phi = @(x) opCS(x,1);
        Phit = @(x) opCS(x,2);
        nu_phi = 1;         % Since opRC is a normalized tight frame.
        tight_phi = 1;      % Phi is a tight frame
        Phi1 = Phi;         % Here no difference between Phi1 and Phi. This
                            % is done for compatibilty reason
        Phit1 = Phit;
    else
        error('Error: Unknown sampling method')

    end



    % Measurements without noise
    y = Phi(Img(:)); 


    % Adding noise to the mesaurement depending of the SNR
    sigma = norm(y) * sqrt( 10^(-SNR/10)/(nb_meas*J) );
    z = sigma * randn(nb_meas*J,1);  % create the noise


    y = y + z;              % Additive white gaussian noise.

    % If the measurement are decorrelated, we solve the problem in a
    % subspace. So we replace the measurement y by P_inv(H) * y
    if decorr==1
        % compute P_inv(H) * y
        y = reshape(y, [], J);
        y = y * H*inv(H'*H);
        y = y(:);
    end


%% Define the 3 different prox for the problem  
    
    %****  Proj L2  *******************************************************
    %
    %   Projection on a B2-ball
    %
    %   || Phi (S * H) - y ||_2 < epsilon 
    
    if decorr==1
        % if measurement are decorrelated, we do not need put H in the
        % projection, so we replace it by the eye matrix (compatibility
        % reasons).
        J2 = I;                 % J2 is defined for compatibility reasons
        H1 = eye(I,J2);         % H1 is defined for compatibility reasons
        epsilon = norm( reshape(z,[],J)* H*inv(H'*H),'fro') ; % typicaly z is not known, but epsilon can be estimated with sigma. 
        nu_H = 1;
        tight_H = 1; 
    else
        % if the measurement are correlated, we need to take accound of H.
        % This will slow down drastically the computation since H is not a
        % tight frame
        H1 = H;                 % H1 is defined for compatibility reasons
        J2 = J;                 % J2 is defined for compatibility reasons
        epsilon = norm(z(:));   % typicaly z is not known, but epsilon can be estimated with sigma. 
        tight_H = 0;
        nu_H = norm(H1)^2;
    end   
    
    % Other parameters for the L2 projection
    paramL2.tight = tight_phi * tight_H;
    paramL2.nu = nu_phi * nu_H;
    paramL2.tol = 1e-5;
    paramL2.maxit = 100;
    paramL2.verbose = 1;
    paramL2.y = y;              % measures
    paramL2.epsilon = epsilon;  %radius of the ball
    paramL2.A = @(x) Phi1( reshape(reshape(x,N,I)*H1',[],1 ));      % operator
    paramL2.At = @(x) reshape( reshape(Phit1(x),N,J2)*H1 ,[],1 );   % adjoint operator
  
    
    f1.prox=@(x,T) fast_proj_b2(x,T, paramL2);
    f1.eval=@(x) 0;            
    
    
    
    %****  Proj simplex  **************************************************
    % The home-made simplex_proj function will achieve two projections in
    % order to satify to constraints
    %
    %     -   (S)_{i,j} > 0  for all i,j    (positivity constraint)
    %
    %     -   Sum_j   S_{i,j}  = 1 for all i
    %
    
    param_simplex.algo = 'duchi';
    
    f2.prox = @(x,T) reshape(simplex_proj(reshape(x,N,I),param_simplex),N*I,1);
    f2.eval = @(x) 0;   
    
    
    
    %****  Prox_TV  *******************************************************
    % 
    % Minimization of the TV norm
    %
    %   argmin_S  Sum_j ||S_{.,j}||_TV
    %
    
    param_TV.maxit = 200;
    param_TV.tol = 10e-5;
    param_TV.verbose = 1;
        
    k=1/3;
    
    f3.prox = @(x,T) reshape(prox_tv(reshape(x,n1,[]),k*T,param_TV),N*I,1);    
    f3.eval = @(x) norm_tv(reshape(x,n1,[]));   % This norm is considered
                                                % as the objective function
                                                
    
%     f3.prox = @(x,T) reshape(prox_TV(reshape(x,n1,[]), T, 'beta',0.249,...
%         'min_rel_obj', 10e-5, 'it_max', 200, ...
%         'it_min', 2, 'verbose', 0),N*I,1); 
    
                                          
    %****** Corresponding sparsifying 2D Wavelet Operator *****************
    %   
    %   You can replace the TV minimization constraint by a Wavelet-L1 constraint   
    %
    %   In that case we minimize
    %
    %   argmin_S  Sum_j || W2D' * S_{.,j} ||_1
    %
    %   Need for sparco toolbox !
    
    % Selection of the wavelet
    wname = 'Daubechies'; % 'haar', 'Daubechies'
    
    DecLev = log2(n1);
   
    W2D = opWavelet(n1,n2,wname,8,DecLev,'min'); % Creation of the wavelet

    opW2D_blk = opBlock_diag_same(I,W2D);

    paramW.verbose=1;
    paramW.maxit=100;
    paramW.At=@(x) opW2D_blk(x,1);
    paramW.A=@(x) opW2D_blk(x,2);
    k=0.03333;
    f4.prox=@(x, T) reshape(prox_l1(reshape(x,[],1),k*T,paramW),N*I,1);  
    f4.eval=@(x) norm(paramW.A(x(:)),1);   




%% solve the problem using the toolbox function   
     
    %Initial point
    x_0 = zeros(N*I,1);
    
    %Create the vector of functions
    if strcmp(method ,'TV')
        F = {f1 f2 f3};     % Use the TV norm
        param.gamma = 1;
        param.tol = 3e-5;
    elseif strcmp(method ,'Wavelet-L1')
        F = {f1 f2 f4};    % Use the wavelet
        param.gamma = 0.3 * max(abs(paramW.A(sources(:))));
        param.tol = 1e-6;
    else 
        error('Error: Unknown regularization method')
    end
    % Parameter for ppxa algorithm
    
    param.lambda = 1;
    param.maxit = 180;
    param.verbose = 2;
    % solve the probleme
    [S_est]= ppxa(x_0,F, param);
    
    %reshape the solution
    S_est = reshape(S_est,N,I);
    
    % image estimation with S
    Img_est = S_est*H';

end

