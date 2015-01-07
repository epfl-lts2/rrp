%RR_LRTV_HSI Low rank TV on hyperspectral image
%                           
%   Reproducible research addendum for compressive source separation
%   ----------------------------------------------------------------
%   
%   THEORY AND METHODS FOR HYPERSPECTRAL IMAGING
%
%   Paper: Mohammad Golbabaee and Pierre Vandergheynst
%   
%   Demonstration matlab file:  Perraudin Nathanael, Mohammad Golbabaee
%
%   EPFL -- February 2012
%   
%   Dependencies
%   ------------
%
%   In order to use this matlab file you need the UNLocXbox toolbox. You
%   can download it on https://lts2research.epfl.ch/unlocbox
%   
%   CS recovery of Hyperpectral images via joint Trace/TV norm minimization
%   -----------------------------------------------------------------------
%
%   We would like to solve this problem: 
%
%   ..   argmin_X   ||X||_* + alpha sum_j ||Xj||_TV
%
%   ..   such that  || A (X) - Y ||_F <= epsilon   Projection on a B2-Ball
%
%   .. math:: \begin{split} \operatorname{arg\,min}_X & \|X\|_* + \alpha \sum_j \|X_j\|_{TV}  \\  \text{ such that: } &\| A ( X ) - Y  \|_F < \epsilon  \text{ (Projection on a B2-Ball) }  \end{split}
%
%   with 
%         - $X$ :  N*J,  HSI 
%         - $Y$ =  $A(X) + Z$: the measurement (forward model)
%         - $A$ :  Sampling operator
%
%
%   Techniques for solving the problem
%   ----------------------------------
%
%   * Dense operator with TV regularization (Correlated measurements)
%   
%
%   Results
%   -------
%
%   .. figure::
%
%      Frames 5 and 6
%
%       
%
%   .. figure::
%
%      Frames 32 and 33
%  
%
%
%   References: golbabaee2012compressed golbabaee2010multichannel
%


%
% Author: Nathanael Perraudin, Mohammad Golbabaee
% Feb 2013, Lausanne, Switzerland

%% ****** Initialization ******
clear all
clc
close all
init_unlocbox


%% ****** Loding renorm(stvcpu-stvgpu)/norm(stvcpu)al data and Refinements ******
fprintf('Loading data...\n')

if 0
    load 'MOFFET_256.mat'

    Img = Img(1:20,1:20,1:40); 
    [n1 , n2, J] = size(Img);
    N = n1*n2;
    Img = reshape(Img,[],J);


else
    load 'Geneva.mat'
    Img= sources*H';
    
    % -- reduce the size for speed computation -- %
    Img=reshape(Img,256,256,171);
     Img= Img(1:40,1:40,1:40);
     Img=reshape(Img,[],size(Img,3));
    % -- ------------------------------------- -- %
    
    n1=sqrt(size(Img,1));
    n2=n1;
    N=n2*n1;
    J=size(Img,2);
end

    

%% ****** Compressive Sampling ******
nb_meas = floor(N/8);
CS_mtx = 1;

%Sampling operator identification.
if CS_mtx == 0
    opCS = opDBD_RC_measur_matrix( N,J , nb_meas);   % "Distributed independent" sampling
    %opCS_multichannel = opBlockDiag (ones(J,1), opCS); % "Uniform" sampling   
    %opCS_multichannel = opCS;
    A = @(x) opCS_multichannel(x,1);
    At = @(x) opCS_multichannel(x,2);
    gamma_phi = 1;   % Since opRC is a normalized tight frame.
    
elseif CS_mtx == 1
    opCS = opRC_measur_matrix( N*J, nb_meas*J);   % Dense CS by random convolution(J. Romberg).
    A = @(x) opCS(x,1);
    At = @(x) opCS(x,2);
    gamma_phi = 1;   % opRC is a normalized tight frame.
end


 
%% ****** Sampling the HSI *******
SNR = inf;

fprintf('Taking CS measurements...\n')
y = A(Img(:));
sigma = sqrt( 10^(-SNR/10)/(nb_meas*J) )*norm(y);
z = sigma * randn(nb_meas*J,1); % Additive white gaussian noise.
y = y + z;


%% ****** Recovery of the HSI from CS measurements, using Lowrank-TV approach. ******
fprintf('Recovery commenced...\n')

if (sigma > 0)
    epsilon = sigma * sqrt(nb_meas*J + 2*sqrt(2*nb_meas*J)) ;  % For noisy data (Fidelity bound)
else
    epsilon =0;
end

% -- Use the unlocbox --

% General parameter
verbose=2;                          % Display
rnk = 10;                           % ???
k = .01*N*J;                        % 0???
alpha=0.1*sqrt(2*rnk/k);            % Weight of the TV norm


% parameter for ppxa
param_ppxa.verbose=verbose;         % Display
param_ppxa.tol = 1e-8;              % Tolerance to stop iterating
param_ppxa.maxit = 100;             % Maximum number of iterations
param_ppxa.lambda=1;                % Stepsize
param_ppxa.gamma=0.1*svds(Img,1);   % Regularization parameter

% Projection on a B2 Ball
param_proj.A=A;                     % Forward operator
param_proj.At=At;                   % Adjoint operator
param_proj.y=y;                     % Measurments
param_proj.verbose=verbose-1;       % Display
param_proj.epsilon=epsilon;         % Radius of the 2 Ball

f1.prox = @(x,T) fast_proj_b2(x, T, param_proj);
f1.eval = @(x) norm(A(x)-y)^2;


% TV proximal operator
param_tv.tol=10e-5;                 % Tolerance to stop iterating
param_tv.maxit=200;                 % Maximum number of iterations
param_tv.verbose=verbose-1;         % Display parameter

f2.prox = @(x,T) reshape(prox_tv(reshape(x, n1, n2,J), T*alpha,param_tv),[],1);
f2.eval = @(x) sum(norm_tv(reshape(x,n1,n2,J)));


% Nuclear norm
param_nuclearnorm.verbose=verbose-1;% Display parameter
        
f3.prox = @(x,T) reshape(prox_nuclearnorm(reshape(x, N, J), T, param_nuclearnorm),[],1);
f3.eval = @(x) sum(svd(reshape(x, N, J)));

% Initial point -- zeros
xin=zeros(N*J,1);


% Solve the problem
Img_est = ppxa(xin,{f1,f2,f3}, param_ppxa);

% reform the cube
Img_est = reshape(Img_est, n1,n2,J);
Img = reshape(Img, n1,n2,J);

% Compute the error
rec_err = norm(Img_est(:) - Img(:))/norm(Img(:));

fprintf('The relative reconstruction error is: %g percent for a compression of %g \n',rec_err*100,N/nb_meas);

%% Display some image

% Selected frames to be displayed
frames=[5,6,32,33];

nfigures=ceil(length(frames)/2);

iframes=1;
for ii=1:nfigures
    figure(ii);
    subplot(221);
    imagesc(Img(:,:,frames(iframes)));
    title(strcat('Original image number: ',num2str(frames(iframes))));
    subplot(222);
    imagesc(Img_est(:,:,frames(iframes)));
    title(strcat('Reconstructed image number: ',num2str(frames(iframes))));
    if iframes==length(frames)
        break;
    end
    subplot(223);
    imagesc(Img(:,:,frames(iframes+1)));
    title(strcat('Original image number: ',num2str(frames(iframes+1))));
    subplot(224);
    imagesc(Img_est(:,:,frames(iframes+1)));
    title(strcat('Reconstructed image number: ',num2str(frames(iframes+1))));
    
    iframes=iframes+2;
end




