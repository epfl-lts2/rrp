%RR_IMAGE_SOURCE_SEPARATION_BIG_DATASET Application of the bloc diagonal method on a big dataset
%                           
%   Reproducible research addendum for compressive source separation
%   ----------------------------------------------------------------
%   
%   THEORY AND METHODS FOR HYPERSPECTRAL IMAGING
%
%   Paper: Mohammad Golbabaee, Simon Arberet, and Pierre Vandergheynst
%   
%   Demonstration matlab file: Perraudin Nathanael, Mohammad Golbabaee
%
%   EPFL -- August 2012
%   
%   Dependencies
%   ------------
%
%   In order to use this matlab file you need the UNLocXbox toolbox. You
%   can download it on https://lts2research.epfl.ch/unlocbox/
%
%   The problem
%   -----------
%  
%   We would like to solve this problem: 
%
%   ..   argmin_S   Sum_j ||S_{.,j}||_TV 
%
%   ..   such that  || Phi ( S * H^t ) - Y ||_F < epsilon   Projection on a B2-Ball
%                                                       
%   ..        and   (S)_{i,j} > 0         for all i,j     (positivity constraint)
%
%   ..        and   Sum_j   S_{i,j}  = 1  for all i
%
%   .. math:: \begin{split} \operatorname{arg\,min}_S & \sum_j \|S_p(.,j)\|_{TV} \\  \text{ such that: }  &\| \Phi ( S \cdot H^T ) - y  \|_F < \epsilon  \text{ (Projection on a B2-Ball) } \\ \text{and } &  S(i,j) > 0  \text{ for all } i,j    \text{ (positivity constraint) } \\ \text{and } &  \sum_j   S(i,j)  = 1   \text{ for all } i  \end{split}
%
%   with 
%         - $S$ :  Sources
%         - $Y$ =  $\Phi(I_m) + Z$: the measurement (forward model)
%         - $Z$ :  noise
%         - $I_m$:  Original 3D image  (64 x 64 x 64) in 'Data.mat'
%         - $H$ :  Mixing matrix
%
%   Uniform sampling
%   ----------------
%
%   For UNIFORM sampling (Block_diag), a decorrelation step can be applied
%   by projecting the measurements onto pinv(H), and solve a easier
%   problem:
%   
%   ..   argmin_S   Sum_j ||Sp_{.,j}||_TV 
%       
%   ..   such that  || Phi *  Sp  - Yp ||_F < epsilon  Projection on a B2-Ball
%                                                   
%   ..        and   (Sp)_{i,j} > 0        for all i,j  (positivity constraint)
%
%   ..        and   Sum_j   Sp_{i,j}  = 1 for all i
%   
%   ..        with  Yp =  Y * (P_inv(H))' 
%   
%   .. math:: \begin{split} \operatorname{arg\,min}_S & \sum_j \|S_p(.,j)\|_{TV}   \\ \text{ such that: } &\| \Phi  S_p  - Y_p \|_F < \epsilon  \text{ (Projection on a B2-Ball) } \\ \text{and } &  S_p(i,j) > 0  \text{ for all } i,j   \text{ (positivity constraint) } \\ \text{and } &  \sum_j   S_p(i,j)  = 1   \text{ for all } i \\ \text{with } & Y_p =   Y H^{+,T} \text{ where } (H^+ \text{ denotes the pseudo inverse of }H) \end{split}
%
%   with 
%         - $Y$ =  $\Phi I_m + Z$: the measurement (forward model)
%
%   This step accelerates the computation and increase accuracy of source recovery. 
%   However, the use of the block diagonale (Block_diag) sampling is required!
%
%   Results
%   -------
%
%   .. figure::
%
%      Sprectral signature of the different sources
%
%       
%
%   .. figure::
%
%      Orignial sources 
%
%       
%
%   .. figure::
%
%      Block_diag operator with TV regularization (Uncorrelated measurements)
%
%       
%
%   .. figure::
%
%      Block_diag operator with Wavelets regularization (Uncorrelated measurements)
%
%       
%
%   References: golbabaee2010multichannel
%



%% Initialization

    clear;
    close all;

    % adding path
    addpath(genpath('./'))
    init_unlocbox()
    global GLOBAL_useGPU;
    verbose = 2;
%% General parameter
 
    sampling_ratios = 1/16;             % compression ratio   only 2^(-p) with  p =1,2,3,... (random convolution RC)

    SNR = inf;                          % SNR (dB)-- inf => no noise
    
    sampling_mecanism = 'Block_diag';   
    
    decorr = 1;                         % decorrelation method (ONLY for Block_diag sampling)
    
    
    %reply = input('Would you like TV or Wavelet-L1 minimization? 1 = TV / 2 = Wavelet-L1 [1]: ', 's');
    %if isempty(reply)
    %    reply = '1';
    %end
    
    %if strcmp(reply,'1')
    %   method = 'TV';                      % 'TV' minimization 
    %elseif strcmp(reply,'2')
    %    method = 'Wavelet-L1';             % 'Wavelet-L1' minimization
    %else
    %    error('Wrong input')
    %end
        
        

%% Loading data for the problem

    load 'Data_big.mat'                 % Synthetic Geneve images 64*64*64
    
    Img = sources*H';
    
    [N , J] = size(Img);                % N dimention of the image, J number of image
    
    n1=sqrt(N);                         % square image x and y dimension
    n2=n1;                             
    
    nb_meas = floor(N*sampling_ratios); % Number of measurements per image

    I = size(H,2);                      % Number of expected sources
    
      
    % display the sources and the mixing parameters
    figure(1)
    %subplot(121)
    plot(H)
    title(sprintf(' Mixing Parameters \n (Spectral signatures) \n Sampling ratio: %g \n SNR: %g',sampling_ratios, SNR));
    xlabel('Spectral band')
    ylabel('Reflectance')
    axis([1,size(H,1),min(0,min(min(H)))*1.1,max(max(H))*1.1])
    %subplot(222)
    
        % plot the result in a nice shape
        sources=sources(:);
        sources1=sources(1:N*I/2);
        sources2=sources(1+N*I/2:end);
        sources1=reshape(sources1,n1,[]);
        sources2=reshape(sources2,n1,[]);
    figure(2)
        imagesc([sources1;sources2]);
        set(gca,'xtick',[])

        colormap hot;
        % Plot separtation lines
        hold on;
        N_lines=I/2-1;

        for ii=1:N_lines;
           plot(ii*n2*ones(2*n1,1)+1,1:2*n1,'k','Linewidth', 3);
        end

        plot(1:3*n1,n1*ones(3*n2,1)+1,'k','Linewidth', 3);
        
    title(sprintf(' Ground Truth \n (Orignial sources) S_1,S_2,..., S_I, I=%i ',I)) 
    drawnow;
    
    

%% Compressed sampling step

    %Sampling operator identification.
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
    paramL2.verbose = verbose -1;;
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
    param_TV.verbose = verbose -1;;
    param_TV.useGPU=GLOBAL_useGPU;
        
    k=0.3333;
    
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
   
    W2D = opWavelet(n1,n2,wname,8,DecLev,'min'); % Creation of the wavele

    opW2D_blk = opBlock_diag_same(I,W2D);

    paramW.verbose= verbose -1;
    paramW.maxit=100;
    paramW.At=@(x) opW2D_blk(x,1);
    paramW.A=@(x) opW2D_blk(x,2);
    k=0.033333;
    f4.prox=@(x, T) reshape(prox_l1(reshape(x,[],1),k*T,paramW),N*I,1);  
    f4.eval=@(x) norm(x(:),1);   




%% solve the problem using the toolbox function   
     
    %Initial point
    x_0 = zeros(N*I,1);
    
    %Create the vector of functions

        F = {f1 f2 f3};     % Use the TV norm
        param.gamma = 1;
        param.tol = 3e-5;
        

    % Parameter for ppxa algorithm
    param.maxit = 200;
    param.lambda = 1;
    param.verbose = verbose;
    % TV
    t=cputime;
    % solve the probleme
    [S_est_tv]= ppxa(x_0,F, param);
    estimation_time_tv=cputime-t
    
    %Create the vector of functions
    F = {f1 f2 f4};    
        param.gamma = 0.3 * max(abs(paramW.A(sources(:))));
        param.tol = 1e-6;
    t=cputime;
    % solve the probleme
    [S_est_l1]= ppxa(x_0,F, param);
    estimation_time_l1=cputime-t
    
    %reshape the solution
    S_est_tv = reshape(S_est_tv,N,I);
    S_est_l1 = reshape(S_est_l1,N,I);
    
    % image estimation with S
    Img_est_tv = S_est_tv*H';
    Img_est_l1 = S_est_l1*H';
    


%% Display the result


    % Evalutaion of the error
    Reconstruction_MSE_tv = norm(Img(:)-Img_est_tv(:))/norm(Img(:))
    Sources_MSE_tv = norm(sources(:)-S_est_tv(:))/norm(S_est_tv(:))
    Reconstruction_MSE_l1 = norm(Img(:)-Img_est_l1(:))/norm(Img(:))
    Sources_MSE_l1 = norm(sources(:)-S_est_l1(:))/norm(S_est_l1(:))

   	figure(3)
        S_est=S_est_tv;
        method='TV';
        % plot the result in a nice shape
        S_est=S_est(:);
        S_est1=S_est(1:N*I/2);
        S_est2=S_est(1+N*I/2:end);
        S_est1=reshape(S_est1,n1,[]);
        S_est2=reshape(S_est2,n1,[]);
        imagesc([S_est1;S_est2]);
        set(gca,'xtick',[])
        colormap hot;
        % Plot separtation lines
        hold on;
        N_lines=I/2-1;

        for ii=1:N_lines;
           plot(ii*n2*ones(2*n1,1)+1,1:2*n1,'k','Linewidth', 3);
        end

        plot(1:3*n1,n1*ones(3*n2,1)+1,'k','Linewidth', 3);
    title(sprintf(strcat(' Recovered sources by:',' ',method,' minimization \n Reconstruction MSE (dB): %g \n Sources MSE (dB): %g \n CPU-Time: %g '),20*log10(Reconstruction_MSE_tv),20*log10(Sources_MSE_tv),estimation_time_tv));
   	
    
    figure(4)
        S_est=S_est_l1;
        method='Wavelet-L1';
        % plot the result in a nice shape
        S_est=S_est(:);
        S_est1=S_est(1:N*I/2);
        S_est2=S_est(1+N*I/2:end);
        S_est1=reshape(S_est1,n1,[]);
        S_est2=reshape(S_est2,n1,[]);
        imagesc([S_est1;S_est2]);
        set(gca,'xtick',[])
        colormap hot;
        % Plot separtation lines
        hold on;
        N_lines=I/2-1;

        for ii=1:N_lines;
           plot(ii*n2*ones(2*n1,1)+1,1:2*n1,'k','Linewidth', 3);
        end

        plot(1:3*n1,n1*ones(3*n2,1)+1,'k','Linewidth', 3);
    title(sprintf(strcat(' Recovered sources by:',' ',method,' minimization \n Reconstruction MSE (dB): %g \n Sources MSE (dB): %g \n CPU-Time: %g '),20*log10(Reconstruction_MSE_l1),20*log10(Sources_MSE_l1),estimation_time_l1));
