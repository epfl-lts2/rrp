%ESTIMATION_PSD Estimation of the power spectrum density
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
%   Figure 1 shows the results of our PSD-estimation algorithm on a
%   $10$-nearest neighbors graph of $20'000$ nodes (sensor type) and only
%   $1$ signal. We compare the estimation using frames of $10$, $20$, $30$,
%   $100$ Gaussian filters. $\sigma$ and $\tau$ are adapted to the number
%   of filters. For this experiment $K_2$ is set to $10$ and the Chebysheff
%   polynomial order is $30$ (Except for $M=100$ where we took $100$). The
%   estimated curves are smoothed versions of the PSD. Since the original
%   PSD is smooth, the estimation is sufficient to construct approximate
%   Wiener filters. Note that with $100$ filters, the windows are very
%   concentrated in the spectral domain and broad in the vertex domain.
%   Thus, we loose the averaging effect of the algorithm resulting in a PSD
%   looking like the Fourier transform of the original signal.
%
%   .. figure::
%
%      Results
%
%      PSD estimation on a graph of $20'000$ nodes with $1$ measurements.
%      Our algorithm is able to successively estimate the PSD of a signal.
%
%   References: perraudin2016stationary


% Author : Nathanael Perraudin
% Date: 6 January 2016



%% Initialization
clear
close all;


N = 20000;
Ndata = 1;
Ndatar = 10;
%% Create a graph
gsp_reset_seed(1);

x = rand(N,3);
G = gsp_nn_graph(x);
G = gsp_estimate_lmax(G);


%% Signal creation

% Take a band limited signal
s0 = @(x) sin(x*4*pi/G.lmax).*(x<G.lmax/4);
s = @(x) (s0(x)).^2;
w = randn(N,Ndata); % White noise
x = gsp_filter(G,s0,w); %Original signal

%% PSD estimation

param.method = 'cheby';
param.order = 30;
param.Nrandom = Ndatar;
t = tic;
param.Nfilt = 10;
psd10 = gsp_psd_estimation(G,x,param);
t1 = toc(t);

t = tic;
param.order = 30;
param.Nfilt = 30;
psd30 = gsp_psd_estimation(G,x,param);
t3 = toc(t);
%%
t = tic;
param.order = 30;
param.Nfilt = 100;
psd100 = gsp_psd_estimation(G,x,param);
t4 = toc(t);
%%
% t = tic;
% param.Nfilt = 500;
% psd500 = gsp_psd_estimation(G,x,param);
% t2 = toc(t);
%% Plot the result
lambda = 0 :0.01: G.lmax;
figure(1)
paramplot.position = [100 100 300 200];
h = plot(lambda,s(lambda),'k',...
    lambda,psd100(lambda),'r:', ...
    lambda,psd30(lambda),'g--',... 
    lambda,psd10(lambda),'b'...   
    );
%     lambda,psd500(lambda),'g',...

%
h(1).LineWidth = 2;
h(2).LineWidth = 2;
h(3).LineWidth = 2;
h(4).LineWidth = 2;
% h(5).LineWidth = 2;
hl = legend('Real PSD',...
     ['Est. M = 100 filters - ',num2str(t4,3), ' s.'],...
     ['Est. M = 30 filters - ',num2str(t3,3), ' s.'],... 
     ['Est. M = 10 filters - ',num2str(t1,3), ' s.']...
      );
%      ['Est. M = 500 filters - ',num2str(t2,3), ' s.'],...

set(hl,'Position', [0.3133 0.5175 0.6483 0.3550])
axis tight
xlabel('\lambda (graph spectral domain)')
title(['PSD estimation using ',num2str(Ndata),' signal'])
gsp_plotfig('psd_estimation',paramplot)



