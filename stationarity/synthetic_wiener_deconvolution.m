%SYNTHETIC_WIENER_DECONVOLUTION Wiener deconvolution experiment on a synthetic dataset
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
%   We start with a de-convolution example on a random geometric graph.
%   This can model an array of sensors distributed in space or simply a
%   mesh. The signal is chosen with a low frequency band-limited PSD. To
%   produce the measurements, the signal is convolved with the heat kernel
%   $e^{-\tau x}$. Additionally, we add some Gaussian noise. The heat
%   kernel is chosen because it simulates a diffusion process. Using
%   de-convolution we aim at recovering the original signal before
%   diffusion. For this experiment, we put ourselves in an ideal case and
%   suppose that  both the PSD  of the input signal and the noise level are
%   known.
%
%   Figure 1 and 2 present the results. We observe that Wiener filtering is
%   able to de-convolve the measurements. The second plot shows the
%   reconstruction errors for three different methods: Tikonov presented in
%   problem, TV and Wiener filtering. Wiener filtering performs clearly
%   much better than the other methods because it has a much better prior
%   assumption.
%
%
%   .. figure::
%
%      Signal and filters for a noise level of $0.16$. 
%
%      The convolution kernel is $e^{-\frac{10 x}{\lambda_{\max}}}$.
%
%   .. figure::
%
%      Evolution of the error with respect of the noise.
%
%      Graph de-convolution on a geometric random graph.
%      
%
%   References: perraudin2016stationary



% Author : Nathanael Perraudin
% Date: 6 January 2016

%% Initialization
clear
close all;
gsp_reset_seed(0)


N = 300; % Number of nodes
sigma = 0.1; % Noise level
verbose = 0; % Verbosity level

%% Create a graph

G = gsp_sensor(N);
G = gsp_compute_fourier_basis(G);
G = gsp_adj2vec(G);


%% Create the problem

% Take a band limited filter
s0 = @(x) sin(x*4*pi/G.lmax).*(x<G.lmax/4);
% PSD
s = @(x) (s0(x)).^2;

% Design a covolution kernel
h = gsp_design_heat(G,10);

% Create the wiener filter
wf = @(x) h(x).*s(x)./((h(x)).^ 2 .*s(x)+ sigma^2 + eps);

%% Signal creation for a single experiment

% Original signal
x = gsp_filter(G,s0,randn(N,1)); 
% Convolve the signal
hx = gsp_filter(G,h,x);
% Add noise 
y = hx +sigma*randn(N,1);

% Apply the wiener filter
xbar = gsp_filter(G,wf,y);


%% Plot a single experiment

figure(1)
subplot(221)
gsp_plot_signal(G,x)
title('Original signal');
caxis([min(x) max(x)]);

subplot(222)
gsp_plot_signal(G,y)
caxis([min(x) max(x)]);
title('Measurements');

subplot(223)
gsp_plot_signal(G,xbar)
caxis([min(x) max(x)]);
title('Recovered signal');

subplot(224)
paramplotf.show_sum = 0;
paramplotf.plot_eigenvalues = 0;
wfdisp = @(x) wf(x)/5;
gplot = {h,s,wfdisp};
gsp_plot_filter(G,gplot,paramplotf)
title('Differrent filters');
legobj = legend('Convolution kernel',...
       'PSD of input signal','Wiener (scalled)');
set(legobj,'Position',[0.6667 0.3025 0.2550 0.1175]);
gsp_plotfig('synthetic_deconvolution')


%% Error with respect of the noise

sigma = 0.01:0.01:0.10;

error_y = zeros(length(sigma),1);
error_tik = zeros(length(sigma),1);
error_tv = zeros(length(sigma),1);
error_wiener = zeros(length(sigma),1);
noise = randn(N,1);
param.verbose = verbose;
for ii = 1: length(sigma)
    % Add noise the noise
    y = hx +sigma(ii)*noise;

    % Create the wiener filter
    wf = @(x) h(x).*s(x)./((h(x)).^ 2 .*s(x)+ sigma(ii)^2 + eps);

    % Apply the filter
    xbar = gsp_filter(G,wf,y);


    sol_tik = gsp_tik_deconvolution_noise(G, y, h, sigma(ii), param);
    sol_tv = gsp_tv_deconvolution_noise(G, y, h, sigma(ii), param);

    % Compute the error
    error_y(ii) = norm(y-x)/norm(x);
    error_tik(ii) = norm(sol_tik-x)/norm(x);
    error_tv(ii) = norm(sol_tv-x)/norm(x);
    error_wiener(ii) = norm(xbar-x)/norm(x);

end

%%
figure(2)
paramplot.position = [100,100,300,220];
noise_level = sigma*sqrt(N)/norm(x);
plot(noise_level,error_tik,noise_level,error_tv,noise_level,error_wiener)
xlabel('Noise level');
ylabel('Relative error');
axis tight;
title('Deconvolution relative error')
legend('Tikonov','TV','Wiener','Location','SouthEast');
gsp_plotfig('synthetic_deconvolution_errors',paramplot)