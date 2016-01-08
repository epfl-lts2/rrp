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
%   USPS dataset. We compute the graph using the first $300$ digits and use
%   $4349$ of the remaining ones to test our algorithm. We use a mask
%   covering $50$ per cent of the pixel and various amount of noise. We
%   then average the result over $4349$ experiments (corresponding to the
%   $4349$ digits) to obtain the curves displayed in Figure 2. For this
%   experiment, we also compare to traditional TV de-noising and Tikonow
%   de-noising. The results presented in Figure 2 show that graph
%   optimization is outperforming classical techniques meaning that the
%   grid is not the optimal graph for the USPS dataset. Moreover, Wiener
%   once again outperforms the other graph-based models.
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
%
%   References: perraudin2016stationary


% Author : Nathanael Perraudin
% Date: 6 January 2016



%% Data handling
close all;
clear
Ns = 300;
verbose = 1;
%% Load the data
[x, y] = load_usps_full();
% Data to learn the kernel
X0 = x(:,1:Ns);
% Data to perform the experiment
X = x(:,Ns:end);

%% Graph creation from the data X0
param.use_flann = 1;
param.k = 20;
param.sigma = 0.2*size(x,2);
G = gsp_nn_graph(X0,param);
G = gsp_compute_fourier_basis(G);
G = gsp_adj2vec(G);


%% Covariance matrices

CovM0 = gsp_stationarity_cov(X0);

CovM = gsp_stationarity_cov(X);
CovMF = G.U'*CovM*G.U;

r = gsp_stationarity_ratio(G, CovM);
fprintf('The stationarity ratio is: %d\n', r);

psd = gsp_experimental_psd(G,CovM0);

%%
sigma = mean(sqrt(sum(X.^2)))*(0.02:0.02:0.2);

error_tik = zeros(size(X,2),length(sigma));
error_tv = zeros(size(X,2),length(sigma));
error_tv_classic = zeros(size(X,2),length(sigma));
error_tik_classic = zeros(size(X,2),length(sigma));
error_wiener = zeros(size(X,2),length(sigma));

G2 = gsp_2dgrid(16,16);
G2 = gsp_compute_fourier_basis(G2);

param.verbose = verbose;

for jj = 1:length(sigma)
    jj
    wf = @(x) psd(x)./( psd(x)+ sigma(jj)^2 + eps);

    wl = @(x) sigma(jj).^2./(psd(x)+eps);
    f = @(x) 1./(wl(x)+1);

    

    parfor ii = 1:size(X,2)
        
        Nsig = ii;

        Mask = rand(G.N,1)>0.5;

        s = X(:,Nsig);

        y = s + sigma(jj) * randn(G.N,1);
        y = Mask.*y;
        
        paramproj = struct;
        paramproj.A = @(x) Mask.*x;
        paramproj.At = @(x) Mask.*x;
        paramproj.y = y;
        paramproj.epsilon = sqrt(sum(Mask(:)))*sigma(jj);
        paramproj.verbose = verbose -1;
        ffid = struct;
        ffid.prox = @(x,T) proj_b2(x,T,paramproj);
        ffid.eval = @(x) eps;
        
%         paramprox = struct;
%         paramprox.verbose = 0;
%         paramprox.y = y;
%         ffid2 = struct;
%         ffid2.grad = @(x) 2*Mask.*(Mask.*x-y);
%         ffid2.eval = @(x) norm(Mask.*(x-y)).^2;
%         ffid2.beta = 2;
%         fwiener = struct;
%         fwiener.prox = @(x,T) gsp_filter_analysis(G,f,x);
%         fwiener.eval = @(x) 0.5*norm(wl(G.e).*gsp_gft(G,x))^2;


        ftik_classic = struct;
        ftik_classic.grad = @(x) 2*G2.L*x;
        ftik_classic.eval = @(x) gsp_norm_tik(G2,x);
        ftik_classic.beta = 2*G2.lmax;
        
%         ftik = struct;
%         ftik.grad = @(x) 2*G.L*x;
%         ftik.eval = @(x) gsp_norm_tik(G,x);
%         ftik.beta = 2*G.lmax;    
        
        paramtv_classic = struct;
        paramtv_classic.verbose = verbose -1;
        ftvclassic = struct;
        ftvclassic.prox = @(x,T) reshape(prox_tv(reshape(x,16,16),T,paramtv_classic),[],1);
        ftvclassic.eval = @(x) norm_tv(reshape(x,16,16));
%         
%         paramtv = struct;
%         paramtv.verbose = 0;
%         ftv = struct;
%         ftv.prox = @(x,T) gsp_prox_tv(x,T,G,paramtv);
%         ftv.eval = @(x) gsp_norm_tv(G,x);

        % sol = solvep(y,{ffid,fwiener});
        % sol2 = solvep(y,{ffid,fwiener2});
        paramsolver = struct;
        paramsolver.verbose = verbose;
        sol_tik_classic = solvep(y,{ffid,ftik_classic},paramsolver);
%         sol_tik = solvep(y,{ffid,ftik},paramsolver);
        sol_tv_classic = solvep(y,{ffid,ftvclassic},paramsolver);
%         sol_tv = solvep(y,{ffid,ftv},paramsolver);
        sol_tik = gsp_tik_inpainting_noise(G, y, Mask, sigma(jj), param);
        sol_tv = gsp_tv_inpainting_noise(G, y, Mask, sigma(jj), param);
        A = @(x) Mask.*x;
        At = @(x) Mask.*x;
        sol_wiener = gsp_wiener_inpainting(G,y, A, At, psd, sigma(jj).^2, param);
%         sol_wiener = solvep(y,{ffid2,fwiener},paramsolver);

        error_tik(ii,jj) = norm(sol_tik - s)/norm(s);
        error_tv(ii,jj) = norm(sol_tv - s)/norm(s);
        error_tik_classic(ii,jj) = norm(sol_tik_classic - s)/norm(s);
        error_tv_classic(ii,jj) = norm(sol_tv_classic - s)/norm(s);
        error_wiener(ii,jj) = norm(sol_wiener - s)/norm(s);



    end

    
end
%% Compute mean error
merr_tik = mean(error_tik,1);
merr_tv = mean(error_tv,1);
merr_tik_classic = mean(error_tik_classic,1);
merr_tv_classic = mean(error_tv_classic,1);
merr_wiener = mean(error_wiener,1);   

save('USPS_eperiment.mat')

%% Plot results

figure(1)
paramplot.position = [100,100,600,220];
subplot(121)
imagesc(abs(G.W))
colorbar
title('Graph weighted adjacency matrix');
subplot(122)
a = -10;
disp = 10*log10(abs(CovMF(1:50,1:50)));
disp(disp<a) = a;
imagesc(disp)
colorbar
imagesc(disp)
colorbar
title('Covariance matrix in Fourier (dB)');
gsp_plotfig('usps_cov',paramplot)

nfac = mean(sqrt(sum(X.^2)));
figure(2)
paramplot.position = [100,100,300,220];
plot(sigma/nfac, merr_tik_classic, ...
    sigma/nfac, merr_tv_classic, ...
    sigma/nfac, merr_tik, ...
    sigma/nfac, merr_tv, ...
    sigma/nfac, merr_wiener);
ylabel('Relative error');
xlabel('Noise level');
axis tight;
title('Inpainting relative error')
legend('Classic Tikonov','Classic TV','Graph Tikonov','Graph TV','Wiener','Location','NorthWest');
gsp_plotfig('usps_inpainting_errors',paramplot)



