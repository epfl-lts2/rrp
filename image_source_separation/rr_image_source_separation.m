%RR_IMAGE_SOURCE_SEPARATION Comparison of different methods
%                           
%   Reproducible research addendum for compressive source separation
%   ----------------------------------------------------------------
%   
%   THEORY AND METHODS FOR HYPERSPECTRAL IMAGING
%
%   Paper: Mohammad Golbabaee, Simon Arberet, and Pierre Vandergheynst
%   
%   Demonstration matlab file:  Perraudin Nathanael, Mohammad Golbabaee
%
%   EPFL -- August 2012
%
%   http://infoscience.epfl.ch/record/180911/files/RRarchive_2.zip
%   
%   Dependencies
%   ------------
%
%   In order to use this matlab file you need the UNLocXbox toolbox. You
%   can download it on http://unlocbox.sourceforge.net
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
%   .. math:: \begin{split} \operatorname{arg\,min}_S & \sum_j \|S_p(.,j)\|_{TV}   \\ \text{ such that: } & \| \Phi ( S \cdot H^T ) - y  \|_F < \epsilon  \text{ Projection on a B2-Ball } \\ \text{and } &  S(i,j) > 0  \text{ for all } i,j    \text{ (positivity constraint) } \\ \text{and } &  \sum_j   S(i,j)  = 1   \text{ for all } i  \end{split}
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
%   ..   such that  || Phi * Sp  - Yp ||_F < epsilon  Projection on a B2-Ball
%                                                   
%   ..        and   (Sp)_{i,j} > 0        for all i,j  (positivity constraint)
%
%   ..        and   Sum_j   Sp_{i,j}  = 1 for all i
%   
%   ..        with  Yp =  Y * (P_inv(H))' 
%   
%   .. math:: \begin{split} \operatorname{arg\,min}_S & \sum_j \|S_p(.,j)\|_{TV}  \\  \text{ such that: } & \| \Phi  S_p  - Y_p \|_F < \epsilon \text{ Projection on a B2-Ball } \\ \text{and } &  S_p(i,j) > 0  \text{ for all } i,j    \text{ (positivity constraint) } \\ \text{and } &  \sum_j   S_p(i,j)  = 1   \text{ for all } i \\ \text{with } & Y_p =   Y H^{+,T}  (\text{ where } H^+ \text{ denotes the pseudo inverse of }H) \end{split}
%
%   with 
%         - $Y$ =  $\Phi I_m + Z$: the measurement (forward model)
%
%   This step accelerates the computation and increase accuracy of source recovery. 
%   However, the use of the block diagonale (Block_diag) sampling is required!
%
%   Techniques for solving the problem
%   ----------------------------------
%
%   We will compare 4 different techniques to solve the problem  
%
%   * Block_diag operator with TV regularization (Decorrelated measurements)
%
%   * Block_diag operator with Wavelets regularization (Decorrelated measurements)
%
%   * Dense operator with TV regularization (Correlated measurements)
%
%   * Dense operator with Wavelets regularization (Correlated measurements)
%   
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
%   .. figure::
%
%      Dense operator with TV regularization (Correlated measurements)
%
%       
%
%   .. figure::
%
%      Dense operator with Wavelets regularization (Correlated measurements)
%
%        
%
%   References: golbabaee2010multichannel
%



%% Initialization

    clear all;
    close all;
    clc;

    % adding path
    addpath(genpath('./'))


%% General parameter
 
    sampling_ratios = 1/4;              % compression ratio   only 2^(-p) with  p =1,2,3,... (random convolution RC)

    SNR = inf;                          % SNR (dB)-- inf => no noise
    


%% Loading data for the problem

    load 'Data.mat'                     % Synthetic Geneve images 64*64*64

    [n1, n2 , J] = size(Img);           % n1, n2 : dimention of the image, J number of image
    
    N = n1*n2;                          % Number of pixels per image
    
    nb_meas = floor(N*sampling_ratios); % Number of measurements per image

    I = size(H,2);                      % Number of expected sources
    
      
    % display the sources and the mixing parameters
    figure(1)
    %set(figure(1),'Units','Normalized','OuterPosition',[0 0 1 1])  
    %subplot(121)
    plot(H)
    title(sprintf(' Mixing Parameters \n (Spectral signatures) \n Sampling ratio: %g \n SNR: %g',sampling_ratios, SNR));
    xlabel('Spectral band')
    ylabel('Reflectance')
        axis([1,size(H,1),min(0,min(min(H)))*1.1,max(max(H))*1.1])
    %subplot(522)
    figure(2)
    imagesc(reshape(sources,n1,[]))
    set(gca,'xtick',[])
    colormap hot;
    % Plot separtation lines
    hold on;
    N_lines=I-1;
    
    for ii=1:N_lines;
       plot(ii*n2*ones(n1,1)+1,1:n1,'k','Linewidth', 3);
    end
    
    title(sprintf(' Ground Truth \n (Orignial sources) S_1,S_2,..., S_I, I=%i ',I)) 
    drawnow;
    
%% Method 1: Block_diag (decorrelated measurements) with TV regularization
 
    sampling_mecanism = 'Block_diag';   % 3 possibility 'Dense','DBD','Block_diag'
    
    decorr = 1;                         % decorrelation method (ONLY for Block_diag sampling)
    
    method = 'TV';                      % 'TV' or 'Wavelet-L1' minimization


    t=cputime;
[ S_est, Img_est ] = general_solver( decorr,sampling_mecanism,method,SNR,Img,H,sources, N,n1,n2,I,J, nb_meas );
    estimation_time=cputime-t
    % Evalutaion of the error
    Reconstruction_MSE = norm(Img(:)-Img_est(:))/norm(Img(:))
    Sources_MSE = norm(sources(:)-S_est(:))/norm(S_est(:))
    
    % reshape the solution in "3 images form"
    recovered_sources = reshape(S_est,n1,n2,[]);

    figure(3)
    %subplot(524)
    imagesc(reshape(recovered_sources,n1,[]))
    set(gca,'xtick',[])
    colormap hot;
    % Plot separtation lines
    hold on;
    N_lines=I-1;
    
    for ii=1:N_lines;
       plot(ii*n2*ones(n1,1)+1,1:n1,'k','Linewidth', 3);
    end
    title(sprintf(' Recovered sources - TV - Block-diag '));
    xlabel(sprintf('Reconstruction MSE (dB): %g  Sources MSE (dB): %g  CPU-Time: %g ',20*log10(Reconstruction_MSE),20*log10(Sources_MSE),estimation_time));
    drawnow;
    
   %% Method 2: Block_diag (decorrelated measurements) with Wavelet regularization
 
    sampling_mecanism = 'Block_diag';   % 3 possibility 'Dense','DBD','Block_diag'
    
    decorr = 1;                         % decorrelation method (ONLY for Block_diag sampling)
    
    method = 'Wavelet-L1';                      % 'TV' or 'Wavelet-L1' minimization


    t=cputime;
[ S_est, Img_est ] = general_solver( decorr,sampling_mecanism,method,SNR,Img,H,sources, N,n1,n2,I,J, nb_meas );
    estimation_time=cputime-t
    % Evalutaion of the error
    Reconstruction_MSE = norm(Img(:)-Img_est(:))/norm(Img(:))
    Sources_MSE = norm(sources(:)-S_est(:))/norm(S_est(:))
    
    % reshape the solution in "3 images form"
    recovered_sources = reshape(S_est,n1,n2,[]);

    figure(4)
    %subplot(526)
    imagesc(reshape(recovered_sources,n1,[]))
    set(gca,'xtick',[])
    colormap hot;
    % Plot separtation lines
    hold on;
    N_lines=I-1;
    
    for ii=1:N_lines;
       plot(ii*n2*ones(n1,1)+1,1:n1,'k','Linewidth', 3);
    end
    title(sprintf(' Recovered sources - Wavelet - Block-diag '));
    xlabel(sprintf('Reconstruction MSE (dB): %g Sources MSE (dB): %g CPU-Time: %g ',20*log10(Reconstruction_MSE),20*log10(Sources_MSE),estimation_time));
    drawnow;
    
    %% Method 3: Dense (correlated measurements) with TV regularization
 
    sampling_mecanism = 'Dense';   % 3 possibility 'Dense','DBD','Block_diag'
    
    decorr = 0;                         % decorrelation method (ONLY for Block_diag sampling)
    
    method = 'TV';                      % 'TV' or 'Wavelet-L1' minimization


    t=cputime;
[ S_est, Img_est ] = general_solver( decorr,sampling_mecanism,method,SNR,Img,H,sources, N,n1,n2,I,J, nb_meas );
    estimation_time=cputime-t
    % Evalutaion of the error
    Reconstruction_MSE = norm(Img(:)-Img_est(:))/norm(Img(:))
    Sources_MSE = norm(sources(:)-S_est(:))/norm(S_est(:))
    
    % reshape the solution in "3 images form"
    recovered_sources = reshape(S_est,n1,n2,[]);

    figure(5)
    %subplot(528)
    imagesc(reshape(recovered_sources,n1,[]))
    set(gca,'xtick',[])
    colormap hot;
    % Plot separtation lines
    hold on;
    N_lines=I-1;
    
    for ii=1:N_lines;
       plot(ii*n2*ones(n1,1)+1,1:n1,'k','Linewidth', 3);
    end
    title(sprintf(' Recovered sources - TV - Dense '));
    xlabel(sprintf('Reconstruction MSE (dB): %g  Sources MSE (dB): %g  CPU-Time: %g ',20*log10(Reconstruction_MSE),20*log10(Sources_MSE),estimation_time));
    drawnow;
    
   %% Method 4: Dense (correlated measurements) with Wavelet regularization
 
    sampling_mecanism = 'Dense';   % 3 possibility 'Dense','DBD','Block_diag'
    
    decorr = 0;                         % decorrelation method (ONLY for Block_diag sampling)
    
    method = 'Wavelet-L1';                      % 'TV' or 'Wavelet-L1' minimization

    t=cputime;
[ S_est, Img_est ] = general_solver( decorr,sampling_mecanism,method,SNR,Img,H,sources, N,n1,n2,I,J, nb_meas );
    estimation_time=cputime-t
    % Evalutaion of the error
    Reconstruction_MSE = norm(Img(:)-Img_est(:))/norm(Img(:))
    Sources_MSE = norm(sources(:)-S_est(:))/norm(S_est(:))
    
    % reshape the solution in "3 images form"
    recovered_sources = reshape(S_est,n1,n2,[]);

    figure(6)
    %subplot(5,2,10)
    imagesc(reshape(recovered_sources,n1,[]))
    set(gca,'xtick',[])
    colormap hot;
    % Plot separtation lines
    hold on;
    N_lines=I-1;
    
    for ii=1:N_lines;
       plot(ii*n2*ones(n1,1)+1,1:n1,'k','Linewidth', 3);
    end
    title(sprintf(' Recovered sources - Wavelet - Dense '));
    xlabel(sprintf('Reconstruction MSE (dB): %g  Sources MSE (dB): %g CPU-Time: %g ',20*log10(Reconstruction_MSE),20*log10(Sources_MSE),estimation_time));
    drawnow;
