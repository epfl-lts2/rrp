%EXPERIMENT_MOLENE_HUMIDITY In-painting on the Molene humidity dataset
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
%   Dataset
%   -------
%
%   The French national meteorological service has published in open access
%   a dataset with hourly weather observations collected during the Month
%   of January 2014 in the region of Brest (France). From these data, we
%   wish to ascertain that our method still performs better than the two
%   other models (TV and Tikonov) on real measurements. The graph is built
%   from the coordinates of the weather stations by connecting all the
%   neighbours in a given radius with a weight function $w(i,j) = e^{-d^2\tau}$ 
%   where $\tau$ is adjusted to obtain a average degree around $3$
%   ($\tau$, however, is not a sensitive parameter). For our experiments,
%   we consider every time step as an independent realization of a WSS
%   process. As sole pre-processing, we remove the mean of the temperature.
%   Thanks to the $744$ time observation, we can estimate the covariance
%   matrix and check wether the process is stationary on the graph.
%
%
%   This experiment
%   ---------------
%
%   The result of the experiment with humidity is displayed in Figures 1
%   and 2. The covariance matrix shows a strong correlation between the
%   different weather stations. Diagonalizing it with the Fourier basis of
%   the graph assesses that the meteorological instances are more or less
%   stationary within the distance graph by highlighting its diagonal
%   characteristic. Moreover this diagonal gives us access to the PSD of
%   the data. In our experiment, we solve an in-painting/de-noising problem
%   with a mask operator covering 50 per cent of measurements and various
%   amount of noise. We then average the result over $744$ experiments
%   (corresponding to the $744$ observations) to obtain the curves
%   displayed in Figure 2. We observe that Wiener optimization performs
%   significantly better when the noise level is high and equivalently well
%   to the two other methods for low noise level.
%
%   .. figure::
%
%      Covariance matrices.
%
%      
%
%   .. figure::
%
%      Recovery errors for different noise levels.
%
%      
%
%
%   References: perraudin2016stationary


% Author : Nathanael Perraudin
% Date: 6 January 2016



%% Initialization
close all
clear
gsp_reset_seed(0)
load('meteo_molene_u.mat');

verbose = 0;


%% Prepare the data
x = info{4};
y = info{3};
z = info{5};

zp = 5*z; % Elevation is more important

% Used coordinates
coords = [x,y,zp];

% Remove the mean to the data (alternatively we could remove 273)
X = value;
X = X - mean(X(:));
%% Graph creation

param.k = 5;
param.type = 'knn';
param.rescale = 0;
param.center = 0;
param.epsilon = sqrt(var(coords(:)))*sqrt(3); 
param.sigma = var(coords(:))*sqrt(3)*0.001;
G = gsp_nn_graph(coords,param);
G.plotting.limits = [-92, 2826, 22313, 24769];
G = gsp_compute_fourier_basis(G);
G = gsp_adj2vec(G);

%% Covariance matrix

CovM = gsp_stationarity_cov(X);
% In the spectral domain
covMF = G.U'*CovM*G.U;
r = gsp_stationarity_ratio(G, CovM);
fprintf('The stationarity ratio is: %d\n', r);
%% Compute an estimation of the PSD

psd = gsp_experimental_psd(G,CovM);


%% Experiment
sigma = mean(sqrt(sum(X.^2)))*(0.02:0.02:0.2);

error_tik = zeros(size(X,2),length(sigma));
error_tv = zeros(size(X,2),length(sigma));
error_wiener = zeros(size(X,2),length(sigma));

param.verbose = verbose;

for jj = 1:length(sigma)
    jj

    

    parfor ii = 1:size(X,2)

        Nsig = ii;

        Mask = rand(G.N,1)>0.5;  %#ok<PFBNS>

        s = X(:,Nsig);  %#ok<PFBNS>

        y = s + sigma(jj) * randn(G.N,1); %#ok<PFBNS>
        y = Mask.*y;
        

        % Solve the problems
         sol_tik = gsp_tik_inpainting_noise(G, y, Mask, sigma(jj), param);
         sol_tv = gsp_tv_inpainting_noise(G, y, Mask, sigma(jj), param);
        A = @(x) Mask.*x;
        At = @(x) Mask.*x;
         sol_wiener = gsp_wiener_inpainting(G,y, A, At, psd, sigma(jj).^2, param);
        error_tik(ii,jj) = norm(sol_tik - s)/norm(s);
        error_tv(ii,jj) = norm(sol_tv - s)/norm(s);
        error_wiener(ii,jj) = norm(sol_wiener - s)/norm(s);



    end

    
end
%% Compute mean errors

merr_tik = mean(error_tik,1);
merr_tv = mean(error_tv,1);
merr_wiener = mean(error_wiener,1);

%% Plot results

figure(1)
paramplot.position = [100,100,600,220];
subplot(121)
imagesc(abs(CovM))
colorbar
title('Covariance matrix');
subplot(122)
disp = 10*log10(abs(covMF));
disp(disp<5) = 5;
imagesc(disp)
imagesc(abs(disp))
colorbar
title('Covariance matrix in Fourier (dB)');
gsp_plotfig('humidity_cov',paramplot)
nfac = mean(sqrt(sum(X.^2)));
figure(2)
paramplot.position = [100,100,300,220];
plot(sigma/nfac,merr_tik,sigma/nfac,merr_tv,sigma/nfac,merr_wiener,...
    'LineWidth',2)
ylabel('Relative error');
xlabel('Noise level');
axis tight;
title('Inpainting relative error')
legend('Tikonov','TV','Wiener','Location','NorthWest');
gsp_plotfig('humidity_inpainting_errors',paramplot)




