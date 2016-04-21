%SYNTHETIC_WIENER_INPAINTING Wiener inpainting experiment on a synthetic dataset
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
%   This experiment
%   ---------------
%
%   We use Wiener optimization to solve an in-painting problem. This time,
%   we suppose that the PSD of the input signal is unknown and we estimate
%   it using $50$ signals. Figures 1, 2 and 3 present quantitative results
%   for the in-painting. Again, we compare three different optimization
%   methods: Tikonov, TV and Wiener. Wiener optimization performs much
%   better than traditional methods because the generated data fits the
%   stationarity model. Moreover, we also observe that the PSD
%   approximation does not affect the result.
%
%   .. figure::
%
%      True VS approximated PSD 
%
%      
%
%   .. figure::
%
%      Wiener filters
%
%      
%
%   .. figure::
%
%      In-painting relative error with respect to number of measurements
%
%       
%
%   References: perraudin2016stationary



% Author : Nathanael Perraudin
% Date: 6 January 2016

%% Initialization
clear
close all;
gsp_reset_seed(0)

N = 400; % Size of the graph
K = 1; % Number of signal to estimate the PSD
K2 = 100; % Number of signal to estimate the covariance matrix empirically
M = 50; % Number of experiment
sigma = 0.05; % Noise level
verbose = 0; % Verbosity level

%% Create a graph
G = gsp_sensor(N);
G = gsp_compute_fourier_basis(G);
G = gsp_adj2vec(G);

%% Create a PSD

% Original filter
g = @(x) sin(3*x).*exp(-x);
% PSD
psd = @(x) (g(x)).^2;


%% Generate data

% filter random signals
w = randn(N,K);
x = gsp_filter_analysis(G,g,w);
w2 = randn(N,K2);
x2 = gsp_filter_analysis(G,g,w2);

% Compute an experimental approximation of the PSD
Cov_exp50 = gsp_stationarity_cov(x2);
% psd50 = gsp_experimental_psd(G,Cov_exp);
param_psd.Nfilt = 50;
psd1 = gsp_psd_estimation(G,x,param_psd);


%% Create the wiener filter

h = @(x) 1;

wf = @(x) h(x).*psd(x)./((h(x)).^ 2 .*psd(x)+ sigma^2 + eps);
wf1 = @(x) h(x).*abs(psd1(x))./((h(x)).^ 2 .*abs(psd1(x))+ sigma^2 + eps);

% Discrete filters for plotting

wl = @(x) sigma.^2./(psd(x)+eps);
f = @(x) 1./(wl(x)+1);

wl1 = @(x) sigma.^2./(psd1(x)+eps);
f1 = @(x) 1./(wl1(x)+1);





%% Solve many problem with respect of the number of measurements

percent = 5:5:95;

error_tik = zeros(length(percent),1);
error_tv = zeros(length(percent),1);
error_wiener = zeros(length(percent),1);
error_wiener1 = zeros(length(percent),1);
error_grm = zeros(length(percent),1);
% Generate M new signals
s = gsp_filter(G,g,randn(N,M));
y2 = s+ sigma *randn(N,M);



ind = randperm(N);

param.verbose = verbose;

for ii = 1:length(percent)
    % Create the problem
    y = y2;
    Mask = zeros(N,1);
    Mask(ind(1:round(percent(ii)/100*N)))=1;
    A = @(x) bsxfun(@times, Mask, x);
    At = @(x) bsxfun(@times, Mask, x);
    y = A(y);



    % Solve the problems
    sol_tik = gsp_tik_inpainting_noise(G, y, Mask, sigma, param);
    sol_tv = gsp_tv_inpainting_noise(G, y, Mask, sigma, param);

    sol_wiener = gsp_wiener_l2(G, y, A, At, psd, sigma^2, param);
    sol_wiener1 = gsp_wiener_l2(G, y, A, At, psd1, sigma^2, param);
    sol_grm = grm_estimator(Cov_exp50,Mask,y,sigma^2);
    error_tik(ii) = norm(sol_tik-s,'fro')/norm(s,'fro');
    error_tv(ii) = norm(sol_tv-s,'fro')/norm(s,'fro');

    error_wiener(ii) = norm(sol_wiener-s,'fro')/norm(s,'fro');
    error_wiener1(ii) = norm(sol_wiener1-s,'fro')/norm(s,'fro');
    error_grm(ii) = norm(sol_grm-s,'fro')/norm(s,'fro');

end
%% Plot the result
figure(1);
paramplot.position = [100 100 300 220];
plot(G.e,psd(G.e),G.e,abs(psd1(G.e)),'--','LineWidth',2);
xlabel('Laplacian eigenvalues')
legend('Real PSD','Approximated PSD');

gsp_plotfig('PSD_approx',paramplot);




figure(2);
paramplot.position = [100 100 300 220];
plot(G.e,wf1(G.e),G.e,wf(G.e),'--','LineWidth',2);
xlabel('Laplacian eigenvalues')
legend('Approximated','Exact');
title('Wiener filters')
gsp_plotfig('approx_wiener_filter',paramplot);

figure(3)
paramplot.position = [100,100,450,315];
plot(percent,error_tik,percent,error_tv,percent,error_wiener,percent,error_wiener1,percent,error_grm,'--','LineWidth',2)
xlabel('Percent of measurements');
ylabel('Relative error');
axis tight;
title('In-painting relative error')
legend('Tikonov','TV','Wiener, exact PSD','Wiener, PSD from 1 meas.','Gaussian MAP from 50 meas.');
gsp_plotfig('synthetic_inpainting_errors',paramplot)


%% Solve 1 optimzation problem


% %Create one problem
% % Generate a new signal
% s = gsp_filter(G,g,randn(N,1));
% % Create a mask
% Mask = rand(N,1)>0.7;
% % Add the noise
% y = Mask.*s+ sigma *randn(N,1);
% % Set the unknown element to min(y) for plotting
% y(logical(1-Mask)) = min(y); 


% % Solve the problems
% % sol_tik = gsp_tik_inpainting_noise(G, y, Mask, sigma);
% % sol_tv = gsp_tv_inpainting_noise(G, y, Mask, sigma);
% 
% A = @(x) Mask.*x;
% At = @(x) Mask.*x;
% sol_wiener = gsp_wiener_l2(G, y, A, At, psd, sigma^2);

% figure(4)
% subplot(221)
% gsp_plot_signal(G,s)
% title('Original signal')
% caxis([min(s),max(s)]);
% subplot(222)
% gsp_plot_signal(G,y)
% title('Measurements')
% caxis([min(s),max(s)]);
% subplot(223)
% gsp_plot_signal(G,sol_wiener)
% title('Recovered signal: Wiener')
% caxis([min(s),max(s)]);
% subplot(224)
% gsp_plot_signal(G,sol_tik)
% title('Recovered signal: Tikonov')
% caxis([min(s),max(s)]);
% gsp_plotfig('inpainting_solution')

