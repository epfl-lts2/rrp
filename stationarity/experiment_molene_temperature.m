%EXPERIMENT_MOLENE_TEMPERATURE In-painting on the Molene temperature dataset
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
%   The result of the experiment with temperature is displayed in Figures 1
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
%   .. figure::
%
%      The temperature of the Island of Brehat.
%
%      
%
%   .. figure::
%
%      An example of the process on the graph (first measure).
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
load('meteo_molene_t.mat');
global SAVE

verbose = 0;
gsp_reset_seed(0)
do_experiment = 0;
%% Prepare the data
x = info{4};
y = info{3};
z = info{5};

zp = 5*z; % Elevation is more important

% Used coordinates
coords = [x,y,zp];

% Remove the mean to the data (alternatively we could remove 273)
X = value;
% X = X - mean(X(:));
mX = mean(X,2);
X = X-repmat(mX,1,size(X,2));

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
CovMF = G.U'*CovM*G.U;
r = gsp_stationarity_ratio(G, CovM);
fprintf('The stationarity ratio is: %d\n', r);
%% Compute an estimation of the PSD

psd = gsp_experimental_psd(G,CovM);
% psd = gsp_psd_estimation(G,X);


%% Experiment

if do_experiment
    percent = 0.2:0.1:0.8;
    sigma = 0.2;
    snr_tik = zeros(size(X,2),length(percent));
    snr_tv = zeros(size(X,2),length(percent));
    snr_wiener = zeros(size(X,2),length(percent));

    param.verbose = verbose;



    for jj = 1:length(percent)
        jj



        parfor ii = 1:size(X,2) 

           Nsig = ii;

            Mask = rand(G.N,1)>percent(jj);  %#ok<PFBNS>

            s = X(:,Nsig);  %#ok<PFBNS>

            y = s + sigma * randn(G.N,1); %#ok<PFBNS>
            y = Mask.*y;


            % Solve the problems
            sol_tik = gsp_tik_inpainting_noise(G, y, Mask, sigma, param);
            sol_tv = gsp_tv_inpainting_noise(G, y, Mask, sigma, param);
    %             A = @(x) Mask.*x;
    %             At = @(x) Mask.*x;
    %             sol_wiener = gsp_wiener_l2(G,y, A, At, psd, sigma(jj).^2, param);

            sol_wiener = gsp_wiener_inpainting(G,y,Mask,psd,sigma^2,param)
            snr_tik(ii,jj) = snr(s,sol_tik);
            snr_tv(ii,jj) = snr(s,sol_tv);
            snr_wiener(ii,jj) = snr(s,sol_wiener);
            snr_in(ii,jj) = snr(s(logical(Mask)),y(logical(Mask)));


        end


    end

    %% Compute mean errors

    msnr_tik = mean(snr_tik,1);
    msnr_tv = mean(snr_tv,1);
    msnr_wiener = mean(snr_wiener,1);
    msnr_in = mean(snr_in(~isnan(snr_in)));
    save('data/molene_t_results.mat','msnr_tik','msnr_tv','msnr_wiener','percent','msnr_in');
else
    load('data/molene_t_results.mat')
end

%% Plot results
paramplot.save = SAVE; 
figure(1)
paramplot.position = [100,100,600,220];
subplot(121)
imagesc(abs(CovM))
colorbar
title('Covariance matrix');
subplot(122)
a = -10;
disp = 10*log10(abs(CovMF));
disp(disp<a) = a;
imagesc(disp)
colorbar
imagesc(disp)
colorbar
title('Covariance matrix in Fourier (dB)');
gsp_plotfig('temperature_cov',paramplot)

nfac = mean(sqrt(sum(X.^2)));
figure(2)
paramplot.position = [100,100,300,220];
plot(percent*100,msnr_tik,percent*100,msnr_tv,percent*100,msnr_wiener,...
    'LineWidth',2)
ylabel('Output SNR (dB)');
xlabel('Percent of missing values');
axis tight;
title('Inpainting error')
legend('Tikhonov','TV','Wiener','Location','Best')
gsp_plotfig('temperature_inpainting_errors',paramplot)

figure(3)
plot(lintimeday,value(1,:)-273)
title(['Temperature over time for ',info{2}(1)]);
axis tight
xlabel('Days')
ylabel('Temperature in degree C')
gsp_plotfig('temperature_time',paramplot)


figure(4)
paramplot.position = [100 100 300 220];
gsp_plot_signal(G , value(:,1)-273);
view(2)
title('Temperatures in degree C')
gsp_plotfig('stations_graph',paramplot)


