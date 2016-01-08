%DEMO_GENERATE_SAMPLES Generation of digit
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
%   Let us focus on the digit $3$. For this experiment, we build a $20$
%   nearest neighbors graph with only $50$ samples. Figure 1 and 2
%   shows the eigenvectors of the Laplacian and of the covariance matrix.
%   Because of stationarity, they are very similar. Moreover, they have
%   $3$-like shape. Using the graph and the PSD, it is also possible
%   generate samples by filtering Gaussian random noise with the following
%   PSD based kernel: $g(\lambda_\ell) = \sqrt{\Gamma_{\ell,\ell}}$. The
%   resulting digits have $3$-like shape confirming the that the class is
%   stationary on the nearest neighbors graph.
%
%   .. figure::
%
%      Laplacian eigenvectors associated to the $16$ smallest non-zero eigenvalues
%
%      
%
%   .. figure::
%
%      Covariance eigenvectors associated with the $16$ highest eigenvalues
%
%      
%
%   .. figure::
%
%      Generated samples by filtering Gaussian random noise on the graph.
%
%      
%
%   .. figure::
%
%      PSD of the data (Note the diagonal shape of the matrix)
%
%      
%
%   References: perraudin2016stationary


% Author : Nathanael Perraudin
% Date: 6 January 2016

%% Initialization
close all
clear

gsp_start;
gsp_reset_seed(1)

%% Load the data
nb = 3;  % Digit
Ns = 50; % Number of sample


% Load the data
[x, y,xx, yy] = load_usps_full();
nx = sqrt(size(x,1));
ny = nx;
X = [x, xx];
Y = [y;yy];
X = X(:,Y==nb );
X = X(:,1:Ns);

%% Create the graph
param.use_full = 1;
param.k = 20;
param.sigma = 0.2*size(X,2); 
G = gsp_nn_graph(X,param);

% Compute the Fourier Basis
G = gsp_compute_fourier_basis(G);

% Compute the covariance matrix
covM = gsp_stationarity_cov(X);

% Compute the covariance matrix in the spectral domain
covMF = G.U'*covM*G.U;

% Compute the covariance eigenvectors
[U,E] = svd(covM);


%% Generate samples
Ng = 25;
w = randn(G.N,Ng);

%gLambda = sqrt(abs(diag(covMF)));
% gs = G.U*(repmat(gLambda,1,size(w,2)).*(G.U'*w));
psd = gsp_experimental_psd(G,covM);
gf = @(x) sqrt(abs(psd(x)));
gs = gsp_filter_analysis(G,gf,w);

% Center and normalize the digits for the plotting
gs = gs-repmat(mean(gs,2),1,Ng);
gs = gs./repmat(std(gs,[],2),1,Ng);



%% Display the results
paramplot.position = [100 100 450 300];
figure()
plot_some_images(U(:,1:16),nx,ny,4,4)
title('Covariance eigenvector')
gsp_plotfig('exp1_cov_eig',paramplot)

figure()
plot_some_images(G.U(:,2:17),nx,ny,4,4)
title('Graph eigenvectors')
gsp_plotfig('exp1_lap_eig',paramplot)

figure
plot_some_images(gs,nx,ny,4,4)
title('Generated samples')
gsp_plotfig('exp1_gen_sample',paramplot)

figure()
imagesc(10*log10(abs(covMF(1:100,1:100))+0.1))
title('PSD matrix in dB');
colorbar
gsp_plotfig('exp1_PSD',paramplot)

%% Compute the stationarity level
r = gsp_stationarity_ratio(G,covM);
fprintf('The stationarity level of the data is: %f\n',r);


