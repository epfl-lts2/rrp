%EXPERIMENT_ORL_INPAINTING In-painting on the ORL dataset
%
%   Authors: Nathanael Perraudin and Pierre Vandergheynst
%
%   Date: January 2016
%
%   Paper: Stationary signal processing on graphs
%
%   Abstract of the paper
%   ---------------------
%   
%   Graphs are a central tool in machine learning and information
%   processing as they allow to conveniently capture the structure of
%   complex datasets. In this context, it is of high importance to develop
%   flexible models of signals defined over graphs or networks. In this
%   paper, we generalize the traditional concept of wide sense stationarity
%   to signals defined over the vertices of arbitrary weighted undirected
%   graphs. We show that stationarity is intimately linked to statistical
%   invariance under a localization operator reminiscent of translation. We
%   prove that stationary signals are characterized by a well-defined Power
%   Spectral Density that can be efficiently estimated even for large
%   graphs. We leverage this new concept to derive Wiener-type estimation
%   procedures of noisy and partially observed signals and illustrate the
%   performance of this new model for denoising and regression.
%
%
%   This experiment
%   ---------------
%
%   For this last experiment, we use the ORL faces dataset. We have a good
%   indication that this dataset is close to stationary since CMUPIE (a
%   smaller faces dataset) is also close to stationary. Each image has
%   $112\times 92=10304$ pixels making it complicated to estimate the
%   covariance matrix and to use a Gaussian MAP estimator. Wiener
%   optimization on the other hand does not necessitate an explicit
%   computation of the covariance matrix. Instead, we estimate the PSD
%   using our aglorithm. A detailed experiment is performed in Figure 3.
%   After adding Gaussian noise to the image, we remove randomly $25\%$ of
%   the pixels. We consider the obtained image as the measurement and we
%   reconstruct the original image using TV, Tikonov and Wiener priors. In
%   Figure 2, we display the reconstruction results for various noise
%   levels. We create the graph with $300$ faces and estimate the PSD with
%   $100$ faces. We test the different algorithms on the $100$ remaining
%   faces.
%
%   .. figure::
%
%      A few images of the dataset
%
%      
%
%   .. figure::
%
%      Recovery errors for different noise levels.
%
%      
%
%   .. figure::
%
%      A single inpainting experiment
%
%      Top left: Original image. Top center: Noisy image (SNR $12.43$ dB).
%      Top right: Measurements $50 \% $ of the noisy image. Bottom left:
%      Reconstruction using Tikhonov prior (SNR $12.12$ dB). Bottom center:
%      Reconstruction using classic TV prior (SNR $13.53$ dB). Bottom
%      right: Reconstruction using Wiener optimization (SNR $14.42$ dB).
%
%
%   References: perraudin2016stationary


% Author : Nathanael Perraudin
% Date: 20 April 2016



%% Data handling
%close all;
clear
gsp_reset_seed(0)
Ng = 300; % Number of samples to construct the graph
Ns = 100; % Number of samples to estimate the PSD
Next = 100; % Number of test samples (All solutions will be computed at once)
verbose = 1; % verbosity
tol =1e-7;  % Tolerance for optimization
maxit = 500; % Maximum number of iterations
perform_simulations = 0;

%% Load the data
[x, y] = load_orl_full();
sx = 112;
sy = 92;
x = x/255;
x = x(:,randperm(400));
% Data to build the graph and learn the kernel
XS = x(:,1:Ng);
% Data to perform the experiment
X = x(:,Ng:(Ng+Next));
mX = mean(X,2);
X = X - repmat(mX,1,size(X,2));
XS = XS - repmat(mX,1,size(XS,2));

% X = X - mean(x(:));
% XS = XS - mean(x(:));
%% Graph creation from the data X0
param.use_flann = 1;
param.k = 10;
%param.sigma = 500;
G = gsp_nn_graph(XS,param);
%G = gsp_2dgrid(112,92);
%parampatch.nnparam = param;
%G = gsp_patch_graph(reshape(XS,sx,sy,Ng),parampatch);
% G = gsp_create_laplacian(G,'combinatorial');
G = gsp_estimate_lmax(G);
% G = gsp_adj2vec(G);
%%
if perform_simulations
    %% Estimate the PSD
    param.order = 100;
    param.Nfilt = 200;
    psd = gsp_psd_estimation(G,XS(:,1:Ns),param);
    %%

%     rel_sigma = (0.05:0.05:0.5);
%     sigma = norm(X,'fro')/sqrt(numel(X))* rel_sigma;
    sigma = 0.05;
    percent = 0.25:0.1:0.75;
    snr_tv_classic = zeros(size(X,2),length(sigma));
    snr_tik = zeros(size(X,2),length(sigma));
    snr_wiener = zeros(size(X,2),length(sigma));    
    snr_tv_classic2 = zeros(size(X,2),length(sigma));
    snr_tik2 = zeros(size(X,2),length(sigma));
    snr_wiener2 = zeros(size(X,2),length(sigma));
    snr_y = zeros(size(X,2),length(sigma));


    param.verbose = verbose;
    % parpool(6)
    for jj = 1:length(percent)



         for ii = 1:size(X,2)
             if verbose
                 fprintf(['Experiment numbers: ',num2str(jj/length(percent)) '   ,   ',num2str(ii/size(X,2)),'\n']);
             end
             s = X(:,ii);
            Mask = rand(G.N,1)>percent(jj);
            A = @(x) bsxfun(@times, Mask, x);
            At = @(x) bsxfun(@times, Mask, x);
            

            y = s + sigma * randn(G.N,size(s,2));
            y = A(y);

            % Classic solution
            paramproj = struct;
            paramproj.A = A;
            paramproj.At = At;
            paramproj.y = y;
            paramproj.epsilon = sqrt(sum(Mask(:)))*sigma;
            paramproj.tight = 1;
            paramproj.verbose = verbose -1;
            ffid = struct;
            ffid.prox = @(x,T) proj_b2(x,T,paramproj);
            ffid.eval = @(x) eps;

            paramtv_classic = struct;
            paramtv_classic.verbose = verbose -1;
            ftvclassic = struct;
            ftvclassic.prox = @(x,T) reshape(prox_tv(reshape(x+mX,sx,sy,size(x,2)),T,paramtv_classic),sx*sy,size(x,2))-mX;
            ftvclassic.eval = @(x) sum(norm_tv(reshape(x+mX,sx,sy,size(x,2))));

            paramsolver = struct;
            paramsolver.verbose = verbose;
            paramsolver.tol = tol;
            paramsolver.maxit = maxit;
            paramsolver.gamma = 0.1;
            sol_tv_classic = solvep(y,{ffid,ftvclassic},paramsolver);

            % Graph solution
            paramsolver = struct;
            paramsolver.verbose = verbose;
            paramsolver.tol = tol;
            paramsolver.order = 40;
            paramsolver.maxit = maxit;
            sol_tik = gsp_tik_inpainting_noise(G, y, Mask, sigma, paramsolver);




            sol_wiener = gsp_wiener_inpainting(G,y,Mask,psd,sigma^2,paramsolver);
            
            snr_mask = @(x,y,Mask) snr(x(logical(1-Mask)),y(logical(1-Mask)));
            snr_tik2(ii,jj) = snr_mask(s,sol_tik,Mask);
            snr_tv_classic2(ii,jj) = snr_mask(s,sol_tv_classic,Mask);
            snr_wiener2(ii,jj) = snr_mask(s,sol_wiener,Mask);
            snr_tik(ii,jj) = snr(s,sol_tik);
            snr_tv_classic(ii,jj) = snr(s,sol_tv_classic);
            snr_wiener(ii,jj) = snr(s,sol_wiener);
            snr_y(ii,jj) = snr(s(logical(Mask)),y(logical(Mask)));
            snr(s,sol_tv_classic)
            snr(s,sol_wiener)

         end


    end
    % Compute mean error
    msnr_tik = mean(snr_tik,1);
    msnr_wiener = mean(snr_wiener,1);   
    msnr_tv_classic = mean(snr_tv_classic,1);       
    msnr_tik2 = mean(snr_tik2,1);
    msnr_wiener2 = mean(snr_wiener2,1);   
    msnr_tv_classic2 = mean(snr_tv_classic2,1);   
    msnr_y = mean(snr_y,1);   

    save('data/ORL_experiment_f.mat','msnr_tik','msnr_wiener','msnr_tv_classic','percent','psd',...
        'msnr_tik2','msnr_wiener2','msnr_tv_classic2','snr_y');
    
else
    load ORL_experiment.mat
end

%% Plot results

paramplot.position = [100,100,300,220];

figure(1)
plot_some_images(X(:,1:6),sx,sy,2,3)
gsp_plotfig('orl_example',paramplot)

%%
paramplot.position = [100,100,300,220];

figure(2)
plot(percent, msnr_tik2, ...
    percent, msnr_tv_classic2, ...
    percent, msnr_wiener2,...
    'LineWidth',2);
ylabel('Ouput SNR (dB)');
xlabel('Percent of missing values');
axis tight;
title('ORL inpainting')
legend('Graph Tikonov','Classic TV','Wiener optimization','Location','Best');
gsp_plotfig('orl_inpainting_errors',paramplot)

%% Perform a simple experiment

if ~perform_simulations
    %% Estimate the PSD
    param.order = 100;
    param.Nfilt = 200;
    psd = gsp_psd_estimation(G,XS(:,1:Ns),param);
end

s = X(:,3);
Mask = rand(G.N,1)>0.5;
A  = @(x) bsxfun(@times, Mask, x);
At = @(x) bsxfun(@times, Mask, x);
sigma_u = 0.05;
y1 = s + sigma_u * randn(G.N,size(s,2));
y = A(y1);

% Classic solution
paramproj = struct;
paramproj.A = A;
paramproj.At = At;
paramproj.y = y;
paramproj.epsilon = sqrt(sum(Mask(:)))*sigma_u;
paramproj.verbose = verbose -1;
paramproj.tight = 1;
ffid = struct;
ffid.prox = @(x,T) proj_b2(x,T,paramproj);
ffid.eval = @(x) eps;

paramtv_classic = struct;
paramtv_classic.verbose = verbose -1;
ftvclassic = struct;
ftvclassic.prox = @(x,T) reshape(prox_tv(reshape(x+mX,sx,sy,size(x,2)),T,paramtv_classic),sx*sy,size(x,2))-mX;
ftvclassic.eval = @(x) sum(norm_tv(reshape(x+mX,sx,sy,size(x,2))));

paramsolver = struct;
paramsolver.verbose = verbose;
paramsolver.gamma = 0.1;
paramsolver.tol = 1e-10;
paramsolver.maxit = maxit;
sol_tv_classic = solvep(y,{ffid,ftvclassic},paramsolver);

% Graph solution
paramsolver = struct;
paramsolver.verbose = verbose;
paramsolver.tol = 1e-10;
paramsolver.order = 40;
paramsolver.maxit = maxit;
sol_tik = gsp_tik_inpainting_noise(G, y, Mask, sigma_u, paramsolver);


sol_wiener = gsp_wiener_inpainting(G,y,Mask,psd,sigma_u^2,paramsolver);

snr(s,y1)
snr(s,sol_tik)
snr(s,sol_tv_classic)

snr(s,sol_wiener)

%%

y(~Mask) = min(y); 
fig = figure();
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.01 0.01], [0.01 0.01]);

subplot(2,3,1);
imagesc(reshape(s+mX,sx,sy))
colormap gray;
axis off;
axis equal
% title('(a)');
subplot(2,3,2);
imagesc(reshape(y1+mX,sx,sy))
colormap gray;
axis off;
axis equal
% title('(b)');
subplot(2,3,3);
imagesc(reshape(y+Mask.*mX,sx,sy))
colormap gray;
axis off;
axis equal
% title('(c)');
subplot(2,3,4);
imagesc(reshape(sol_tik+mX,sx,sy))
colormap gray;
axis off;
axis equal
% title('(d)')
subplot(2,3,5);
imagesc(reshape(sol_tv_classic+mX,sx,sy))
colormap gray;
axis off;
axis equal
% title('(e)')
subplot(2,3,6);
imagesc(reshape(sol_wiener+mX,sx,sy))
colormap gray;
axis off;
axis equal
% title('(f)')
% tightfig(fig);
clear subplot
paramplot.position = [100,100,600,440];
gsp_plotfig('orl_single',paramplot)


