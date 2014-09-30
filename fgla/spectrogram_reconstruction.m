function [ f,Sp,info_reconstruct] = spectrogram_reconstruction( S,G,Gt,...
                                     method,alpha, param)
%SPECTROGRAM_RECONSTRUCTION Reconstruct a signal from a spectrogram
%   Usage: sp = spectrogram_reconstruction(S,A,G,Gt, method, param);
%          [ sp,Sp,info_reconstruct] = spectrogram_reconstruction(S,A,G,... 
%                                       Gt, method, param);
%
%   Input parameters:
%         S     : Spectrogram
%         G     : Analysis operator
%         Gt    : Dual sythesis operator
%         method: Reconstruction method
%         param : Structure of optional parameter
%   
%   Output parameters:
%         f     : reconstructed signal.
%         Sp    : Projected reconstructed signal
%         info_reconstruct: information for the reconstruction
%
%   This function reconstruct a signal $f$ from its spectrogram $S$. 
%   
%   ..        S = |G f|
%
%   .. math   S = |G f|    
%
%   If $S$ has a phase, it will be used as starting phase for the
%   algorithm.
%
%   *method* define the algorithm used to solve the problem. 
%   
%   * 'GLA'             : Griffin-lim algorithm
%   * 'FGLA'            : Fast Griffin-Lim algorithm
%
%   *param* is an optional structure containing optional parameter:
%
%   * *param.maxit*     : Maximum number of iteration (default 100)
%   * *param.tol*       : Tolerance for convergence (default 0.001)
%   * *param.verbose*   : Display parameter (0 no log, 1 main steps, 2
%                         display figure) (default 0)
%   * *alpha*     : Parameter for Forward-PBL and Backwar-PBL.
%                         This must be between 0 and 1 (Default 0.99)
%   * *param.zero_phase*: Starting the algorithm with a zero phase 
%                         (default 0)
%   * *param.rand_phase*: Starting the algorithm with a random phase 
%                         (default 0)
%
%   *info_reconstruct* is a structure containing convergence information.
%
%   * *info_reconstruct.stopping_criterion* : stopping criterion
%   * *info_reconstruct.n_iter*             : final number of iteration
%   * *info_reconstruct.n_cputime*          : cpu time of computation
%                                             This algo should be twice
%                                             faster.
%   * *info_reconstruct.tol*                : final relative difference
%                                             between two iterations
%   * *info_reconstruct.ssnr*               : ssnr each iteration
%




% Optional input arguments
if nargin<6, param=struct; end
if ~isfield(param, 'maxit'), param.maxit = 100; end
if ~isfield(param, 'tol'), param.tol = 0.001; end
if ~isfield(param, 'verbose'), param.verbose = 0; end
if ~isfield(param, 'zero_phase'), param.zero_phase = 0; end
if ~isfield(param, 'random_phase'), param.random_phase = 0; end
if nargin<5, alpha=0.99; end

if nargin<4
    method = 'FGLA';
end


% the initial phase is not known
if param.zero_phase
   S=abs(S); 
end


% the initial phase is random
if param.random_phase
   S=abs(S).*exp(1i*2*pi*rand(size(S))); 
end

% Define the projections
objS=abs(S);    
proj_c1=@(x) G(Gt(x));
proj_c2=@(x) objS.*exp(1i*angle(x));
    


% Start the counter. However, the code is not optimized due to the
% computation of the error every iteration. It could be twice faster
tic;


% Initialisation
Sp=S;

% Structure to store informations    
info_reconstruct=struct;
info_reconstruct.ssnr=zeros(param.maxit,1);
ii=0;
% Griffin method
if strcmp(method,'GLA')
    if param.verbose>0
        fprintf('   * Selected method: Griffin-Lim algorithm \n');
    end

    % Loop
    while 1
        % Update rule   
        ii=ii+1;                    % update current iteration number
        Sp_old=Sp;                  % save the previous result
        Sp=proj_c1(proj_c2(Sp));    % perform the algorithm
        
        % test convergence
        [stop,info_reconstruct]=test_convergence(Sp,Sp_old,objS,ii, ...
                        info_reconstruct,proj_c1,param);
                    
        if stop
            break;   
        end  

    end
    
    


    

    
elseif strcmp(method,'FGLA')
    fprintf('   * Selected method: fast Griffin-Lim algorithm\n');
    
    % Initialisation 
    Tp_old=proj_c1(proj_c2(Sp));
    
    while 1
        % Update rule
        ii=ii+1;                    % update current iteration number
        Sp_old=Sp;                  % save the previous result
        Tp=proj_c1(proj_c2(Sp));    % perform the algorithm
        Sp=Tp+alpha*(Tp-Tp_old);
        Tp_old=Tp;                  % update of Tp_old        
        
        % test convergence
        [stop,info_reconstruct]=test_convergence(Sp,Sp_old,objS,ii, ...
                        info_reconstruct,proj_c1,param);
                    
        if stop
            break;   
        end  

    end 


% elseif strcmp(method,'Backward-PBL')
%     fprintf('   * Selected method: Backward-PBL\n');
%     
%     % Initialisation 
%     Tp_old=proj_c1(proj_c2(Sp));
%     
%     while 1
%         % Update rule
%         Sp_old=Sp;                  % save the previous result
%         ii=ii+1;                    % update current iteration number
%         Tp=proj_c1(proj_c2(Sp));    % perform the algorithm
%         Sp=Tp+alpha*(Tp-Tp_old);
%         Tp_old=Tp;                  % update of Tp_old
%         
%         % test convergence
%         [stop,info_reconstruct]=test_convergence(Sp,Sp_old,objS,ii, ...
%                         info_reconstruct,proj_c1,param);
%                     
%         if stop
%             break;   
%         end  
% 
%     end 


else 
   error('Unknown reconstruction method!');
end

f=Gt(Sp);
Sp=G(f);


end



function  [stop,info_reconstruct]=test_convergence(Sp,Sp_old,objS,ii,info_reconstruct,proj_c1,param)
    % computation for iteration ii
    tol_obs=sum(abs(Sp(:)-Sp_old(:)))/sum(abs(Sp(:)));
    ssnrm=ssnr(objS,proj_c1(Sp));
    info_reconstruct.ssnr(ii)=ssnrm;
    
    % display the result
    if param.verbose>0
        fprintf('       Iteration %i',ii);
        fprintf('     > eps: %g      ssnr: %g \n',tol_obs,ssnrm);
    end

    % test if the algorithm has converged
    stop=0;
    if tol_obs<param.tol
        if param.verbose>0
            fprintf('   * Stopping criterion: epsilon with: tol_obs = %g \n',tol_obs);
        end
        info_reconstruct.stopping_criterion='epsilon';
        stop=1;
    elseif ii==param.maxit
        if param.verbose>0
            fprintf('   * Stopping criterion: max iteration with: tol_obs = %g \n',tol_obs);
        end
        stop=1;
    end
    
    % if we have converged
    if stop
        info_reconstruct.stopping_criterion='max_it';
        info_reconstruct.n_iter=ii;
        info_reconstruct.cputime=toc;
        info_reconstruct.tol=tol_obs; 
    end
    
    
    % Plot evolution on image
    if param.verbose>1
       figure(25)
       subplot(221)
       imagesc(abs(Sp));
       axis('xy');
       title('Modulus')
       subplot(222)
       imagesc(angle(Sp));
       title('Phase')
       axis('xy');
       subplot(223)
       imagesc(abs(Sp)-abs(Sp_old));
       title('Difference Sp Sp-old');
       axis('xy');
%        colorbar;
       drawnow;
    end
end


