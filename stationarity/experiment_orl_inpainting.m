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
%      Top left: Original image. Top center: Noisy image. Top right:
%      Measurements $75 \% $ of the noisy image. Bottom left:
%      Reconstruction using Tikonov prior (relative error $21.3\%$). Bottom
%      center: Reconstruction using classic TV prior (relative error
%      $19.7\%$). Bottom right: Reconstruction using Wiener optimization
%      (relative error $18.2\%$).
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
verbose = 0; % verbosity
tol =1e-7;  % Tolerance for optimization
maxit = 1000; % Maximum number of iterations
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
X = X - mean(x(:));
XS = XS - mean(x(:));
%% Graph creation from the data X0
param.use_flann = 1;
param.k = 10;
%G = gsp_nn_graph(XS,param);
%G = gsp_2dgrid(112,92);
parampatch.nnparam = param;
 G = gsp_patch_graph(reshape(XS,sx,sy,Ng),parampatch);
G = gsp_estimate_lmax(G);
% G = gsp_adj2vec(G);


%% Estimate the PSD

param.order = 100;
param.Nfilt = 200;
psd = gsp_psd_estimation(G,XS(:,1:Ns),param);
%%

rel_sigma = (0.05:0.05:0.5);
sigma = norm(X,'fro')/sqrt(numel(X))* rel_sigma;
% sigma = mean(sqrt(sum(X.^2)))*(0.02:0.04:0.2);

if perform_simulations
error_tv_classic = zeros(size(X,2),length(sigma));
error_tik = zeros(size(X,2),length(sigma));
error_wiener = zeros(size(X,2),length(sigma));


param.verbose = verbose;
% parpool(6)
for jj = 1:length(sigma)
    


     parfor ii = 1:size(X,2)
         if verbose
             fprintf(['Experiment number',num2str(jj/length(sigma)) '   ,   ',num2str(ii/size(X,2)),'\n']);
         end
         s = X(:,ii);
        Mask = rand(G.N,1)>0.25;
        A = @(x) bsxfun(@times, Mask, x);
        At = @(x) bsxfun(@times, Mask, x);


        y = s + sigma(jj) * randn(G.N,size(s,2));
        y = A(y);
        
        % Classic solution
        paramproj = struct;
        paramproj.A = A;
        paramproj.At = At;
        paramproj.y = y;
        paramproj.epsilon = sqrt(sum(Mask(:)))*sigma(jj);
        paramproj.verbose = verbose -1;
        ffid = struct;
        ffid.prox = @(x,T) proj_b2(x,T,paramproj);
        ffid.eval = @(x) eps;
         
        paramtv_classic = struct;
        paramtv_classic.verbose = verbose -1;
        ftvclassic = struct;
        ftvclassic.prox = @(x,T) reshape(prox_tv(reshape(x,sx,sy,size(x,2)),T,paramtv_classic),sx*sy,size(x,2));
        ftvclassic.eval = @(x) sum(norm_tv(reshape(x,sx,sy,size(x,2))));
        
        paramsolver = struct;
        paramsolver.verbose = verbose;
        paramsolver.gamma = 0.1;
        paramsolver.tol = tol;
        paramsolver.maxit = maxit;
        sol_tv_classic = solvep(y,{ffid,ftvclassic},paramsolver);
        
        % Graph solution
        paramsolver = struct;
        paramsolver.verbose = verbose;
        paramsolver.tol = tol;
        paramsolver.order = 100;
      paramsolver.maxit = maxit;
        sol_tik = gsp_tik_inpainting_noise(G, y, Mask, sigma(jj), paramsolver);
        
        
        sol_wiener = gsp_wiener_l2(G,y, A, At, psd, sigma(jj).^2, paramsolver);
        error_tik(ii,jj) = norm(sol_tik - s,'fro')/norm(s,'fro');
        error_wiener(ii,jj) = norm(sol_wiener - s,'fro')/norm(s,'fro');
        error_tv_classic(ii,jj) = norm(sol_tv_classic - s,'fro')/norm(s,'fro');


     end

    
end
% Compute mean error
merr_tik = mean(error_tik,1);
merr_wiener = mean(error_wiener,1);   
merr_tv_classic = mean(error_tv_classic,1);   

save('ORL_eperiment.mat','merr_tik','merr_wiener','merr_tv_classic','rel_sigma');
else
    load('ORL_experiment.mat');
end

%% Plot results


paramplot.position = [100,100,300,220];

figure(1)
plot_some_images(X(:,1:6),sx,sy,2,3)
gsp_plotfig('orl_example',paramplot)

%%
paramplot.position = [100,100,300,220];

figure(2)
plot(rel_sigma, merr_tik, ...
    rel_sigma, merr_tv_classic, ...
    rel_sigma, merr_wiener,...
    'LineWidth',2);
ylabel('Relative error');
xlabel('Noise level');
axis tight;
title('Inpainting relative error')
legend('Graph Tikonov','Classic TV','Wiener optimization','Location','SouthEast');
gsp_plotfig('orl_inpainting_errors',paramplot)

%% Perform a simple experiment

s = X(:,5);
Mask = rand(G.N,1)>0.25;
A = @(x) bsxfun(@times, Mask, x);
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
ffid = struct;
ffid.prox = @(x,T) proj_b2(x,T,paramproj);
ffid.eval = @(x) eps;

paramtv_classic = struct;
paramtv_classic.verbose = verbose -1;
ftvclassic = struct;
ftvclassic.prox = @(x,T) reshape(prox_tv(reshape(x,sx,sy,size(x,2)),T,paramtv_classic),sx*sy,size(x,2));
ftvclassic.eval = @(x) sum(norm_tv(reshape(x,sx,sy,size(x,2))));

paramsolver = struct;
paramsolver.verbose = verbose;
paramsolver.gamma = 0.1;
paramsolver.tol = 1e-10;
paramsolver.maxit = 1000;
sol_tv_classic = solvep(y,{ffid,ftvclassic},paramsolver);

% Graph solution
paramsolver = struct;
paramsolver.verbose = verbose;
paramsolver.tol = 1e-10;
paramsolver.order = 100;
paramsolver.maxit = 1000;
sol_tik = gsp_tik_inpainting_noise(G, y, Mask, sigma_u, paramsolver);


sol_wiener = gsp_wiener_l2(G,y, A, At, psd, sigma_u.^2, paramsolver);
norm(sol_tik - s,'fro')/norm(s,'fro')
norm(sol_tv_classic - s,'fro')/norm(s,'fro')

norm(sol_wiener - s,'fro')/norm(s,'fro')

%%

increase_v =0.06;
y(~Mask) = min(y); 
fig = figure();

h = subplot(231);
p = get(h, 'pos');
p(3:4) = p(3:4) +increase_v;
set(h,'pos',p);
imagesc(reshape(s,sx,sy))
colormap gray;
axis off;
axis equal
% title('(a)');
h = subplot(232);
p = get(h, 'pos');
p(3:4) = p(3:4) +increase_v;
set(h,'pos',p);
imagesc(reshape(y1,sx,sy))
colormap gray;
axis off;
axis equal
% title('(b)');
h = subplot(233);
p = get(h, 'pos');
p(3:4) = p(3:4) +increase_v;
set(h,'pos',p);
imagesc(reshape(y,sx,sy))
colormap gray;
axis off;
axis equal
% title('(c)');
h = subplot(234);
p = get(h, 'pos');
p(3:4) = p(3:4) +increase_v;
set(h,'pos',p);
imagesc(reshape(sol_tik,sx,sy))
colormap gray;
axis off;
axis equal
% title('(d)')
h = subplot(235);
p = get(h, 'pos');
p(3:4) = p(3:4) +increase_v;
set(h,'pos',p);
imagesc(reshape(sol_tv_classic,sx,sy))
colormap gray;
axis off;
axis equal
% title('(e)')
h = subplot(236);
p = get(h, 'pos');
p(3:4) = p(3:4) +increase_v;
set(h,'pos',p);
imagesc(reshape(sol_wiener,sx,sy))
colormap gray;
axis off;
axis equal
% title('(f)')
tightfig(fig);

paramplot.position = [100,100,600,440];
gsp_plotfig('orl_single',paramplot)


