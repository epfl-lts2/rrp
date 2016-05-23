%EXPERIMENT_USPS_INPAINTING In-painting on the USPS dataset
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
%   We perform the same kind of in-painting/de-noising experiment with the
%   USPS dataset. For our experiments, we consider every digit as an
%   independent realization of a GWSS process. As sole pre-processing, we
%   remove the mean (over pixels and digits). It is also possible to remove
%   the mean of each pixel separately. It might increase the stationarity
%   level of the data. In this contribution, we choose not to perform this
%   pre-processing as we consider the raw data stationary. We create the
%   graph using patches of pixels of size $5 \times 5$. The pixels' patches
%   help because we have only a few digits available. When the size of the
%   data increases, a nearest neighbor graph performs even better. We
%   estimate the PSD using only the first $20$ digits and we use $500$ of
%   the remaining ones to test our algorithm. We use a mask covering $50$
%   per cent of the pixel and various amount of noise. We then average the
%   result over $500$ experiments (corresponding to the $500$ digits) to
%   obtain the curves displayed in Figure 2. For this experiment, we also
%   compare to traditional TV de-noising and Tikonov de-noising.
%   Additionally we compute the classical MAP estimator based on the
%   empirical covariance matrix for the solution. The results presented in
%   Figure 2 show that graph optimization is outperforming classical
%   techniques meaning that the grid is not the optimal graph for the USPS
%   dataset. Wiener once again outperforms the other graph-based models.
%   Moreover, this experiment shows that our PSD estimation is robust when
%   the number of signals is small. In other words, using the graph allows
%   us for a much better covariance estimation than a simple empirical
%   average. When the number of measurements increases, the MAP estimator
%   improves in performance and eventually outperforms Wiener because the
%   data is close to stationary on the graph.
%
%
%   .. figure::
%
%      NN graph analysis.
%
%      Left Weights matrix of the $20$ nearest neighbor graph (The diagonal
%      shape indicate the grid base topology of the graph). Right: PSD
%      matrix for the first $50$ graph frequencies.
%
%   .. figure::
%
%      Recovery errors for different noise levels.
%
%      Methods using the nearest neighbors graph performs better.
%
%   .. figure::
%
%      Different PSDs
%
%      Compared to $\frac{1}{x}$, the approximation is a smoothed version
%      of the experimental PSD.
%
%   .. figure::
%
%      Some digits of the USPS dataset.
%
%      
%
%   References: perraudin2016stationary


% Author : Nathanael Perraudin
% Date: 6 January 2016



%% Data handling
close all;
clear
gsp_reset_seed
Ng = 20; % Number of samples to construct the graph
Ns = 20; % Number of samples to estimate the PSD
Next = 500; % Number of test samples
verbose = 0; % verbosity
perform_simulations = 0;
%% Load the data
[x, y] = load_usps_full();
% Data to learn the kernel
X0 = x(:,1:Ns);
% Data to learn the kernel
XG = x(:,1:Ng);
% Data to perform the experiment
X = x(:,Ns:(Ns+Next));
X = X - mean(X(:));
%% Graph creation from the data X0
% param.use_flann = 0;
% param.k = 20;
% param.sigma = 0.2*size(XG,2);
% G = gsp_nn_graph(XG,param);
% G = gsp_2dgrid(16,16);
G = gsp_patch_graph(reshape(XG,16,16,Ng));
G = gsp_compute_fourier_basis(G);
G = gsp_adj2vec(G);

X0 = X0 - mean(X(:));
XG = XG - mean(X(:));
%% Covariance matrices

CovM0 = gsp_stationarity_cov(X0);


CovM = gsp_stationarity_cov(X);
CovMF = G.U'*CovM*G.U;

r = gsp_stationarity_ratio(G, CovM);
fprintf('The stationarity ratio is: %d\n', r);

psd_t = gsp_experimental_psd(G,CovM);
psd = gsp_psd_estimation(G,X0);
%%
rel_sigma = (0.05:0.05:0.5);
sigma = norm(X,'fro')/sqrt(numel(X))*rel_sigma;
% sigma = mean(sqrt(sum(X.^2)))*(0.02:0.04:0.2);

if perform_simulations
    error_tik = zeros(size(X,2),length(sigma));
    error_tv = zeros(size(X,2),length(sigma));
    error_tv_classic = zeros(size(X,2),length(sigma));
    error_tik_classic = zeros(size(X,2),length(sigma));
    error_wiener = zeros(size(X,2),length(sigma));
    error_grm = zeros(size(X,2),length(sigma));

    G2 = gsp_2dgrid(16,16);
    G2 = gsp_compute_fourier_basis(G2);

    param.verbose = verbose;

    for jj = 1:length(sigma)
        jj
    %     wf = @(x) psd(x)./( psd(x)+ sigma(jj)^2 + eps);
    % 
    %     wl = @(x) sigma(jj).^2./(psd(x)+eps);
    %     f = @(x) 1./(wl(x)+1);



        parfor ii = 1:size(X,2)

            Nsig = ii;

            Mask = rand(G.N,1)>0.5;

            s = X(:,Nsig);

            y = s + sigma(jj) * randn(G.N,1);
            y = Mask.*y;

            % Classic solution
            paramproj = struct;
            paramproj.A = @(x) Mask.*x;
            paramproj.At = @(x) Mask.*x;
            paramproj.y = y;
            paramproj.epsilon = sqrt(sum(Mask(:)))*sigma(jj);
            paramproj.verbose = verbose -1;
            ffid = struct;
            ffid.prox = @(x,T) proj_b2(x,T,paramproj);
            ffid.eval = @(x) eps;


            ftik_classic = struct;
            ftik_classic.grad = @(x) 2*G2.L*x;
            ftik_classic.eval = @(x) gsp_norm_tik(G2,x);
            ftik_classic.beta = 2*G2.lmax;


            paramtv_classic = struct;
            paramtv_classic.verbose = verbose -1;
            ftvclassic = struct;
            ftvclassic.prox = @(x,T) reshape(prox_tv(reshape(x,16,16),T,paramtv_classic),[],1);
            ftvclassic.eval = @(x) norm_tv(reshape(x,16,16));


            paramsolver = struct;
            paramsolver.verbose = verbose;
            sol_tik_classic = solvep(y,{ffid,ftik_classic},paramsolver);
            paramsolver.gamma = 0.1;
            sol_tv_classic = solvep(y,{ffid,ftvclassic},paramsolver);

            % Graph solution
            sol_tik = gsp_tik_inpainting_noise(G, y, Mask, sigma(jj), param);
            sol_tv = gsp_tv_inpainting_noise(G, y, Mask, sigma(jj), param);
%             A = @(x) Mask.*x;
%             At = @(x) Mask.*x;
%             sol_wiener = gsp_wiener_l2(G,y, A, At, psd, sigma(jj).^2, param);

            sol_wiener = gsp_wiener_inpainting(G,y,Mask,psd,sigma(jj)^2,param)
            sol_grm = grm_estimator(CovM0,Mask,y,sigma(jj).^2);
            error_tik(ii,jj) = norm(sol_tik - s)/norm(s);
            error_tv(ii,jj) = norm(sol_tv - s)/norm(s);
            error_tik_classic(ii,jj) = norm(sol_tik_classic - s)/norm(s);
            error_tv_classic(ii,jj) = norm(sol_tv_classic - s)/norm(s);
            error_wiener(ii,jj) = norm(sol_wiener - s)/norm(s);

            error_grm(ii,jj) = norm(sol_grm - s)/norm(s);


        end


    end
    % Compute mean error
    merr_tik = mean(error_tik,1);
    merr_tv = mean(error_tv,1);
    merr_tik_classic = mean(error_tik_classic,1);
    merr_tv_classic = mean(error_tv_classic,1);
    merr_wiener = mean(error_wiener,1);   
    merr_grm = mean(error_grm,1);   

    save('USPS_experiment.mat', 'merr_tik', 'merr_tv', 'merr_tik_classic',...
        'merr_tv_classic', 'merr_wiener', 'merr_grm');
else
    load USPS_experiment.mat
end
    
    
%% Plot results

figure(1)
paramplot.position = [100,100,600,220];
subplot(121)
imagesc(abs(G.W))
colorbar
title('Graph weighted adjacency matrix');
subplot(122)
a = -10;
disp = 20*log10(abs(CovMF(1:50,1:50)));
disp(disp<a) = a;
imagesc(disp)
colorbar
imagesc(disp)
colorbar
title('Covariance matrix in Fourier (dB)');
gsp_plotfig('usps_cov',paramplot)
%%
figure(2)
paramplot.position = [100,100,450,330];
plot(rel_sigma, merr_tik_classic, ...
    rel_sigma, merr_tv_classic, ...
    rel_sigma, merr_tik, ...
    rel_sigma, merr_tv, ...
    rel_sigma, merr_wiener,...
    rel_sigma, merr_grm,...
    'LineWidth',2);
ylabel('Relative error');
xlabel('Noise level');
axis tight;
title('In-painting relative error')
legend('Classic Tikonov','Classic TV','Graph Tikonov','Graph TV','Wiener','Gaussian MAP','Location','SouthEast');
gsp_plotfig('usps_inpainting_errors',paramplot)


%% Plot the PSDs
figure(3)
paramplot.position = [100,100,300,220];
plot(G.e,psd_t(G.e),G.e,psd(G.e), G.e, 1./G.e,'LineWidth',2)
axis([0 G.lmax/2 0 15])
h =legend('Experimental (All signals)','Approximation using 20 signals','1/x');
title('PSD')
set(h,'Position', [0.2133 0.6511 0.8200 0.2136])
gsp_plotfig('usps_psd',paramplot)

%% Plot some digits
figure(4)
paramplot.position = [100,100,260,220];
img = flipud(rot90(reshape(X(:,1:16),16,16,16)));
plot_some_images(reshape(img,256,16),16,16,4,4)
title('USPS digits')
gsp_plotfig('usps_digits',paramplot)

