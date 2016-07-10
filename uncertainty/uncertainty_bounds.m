%UNCERTAINTY_BOUNDS Compute the global uncertainty bound of some graphs
%
%   This package contains the code to reproduce all the figures of the
%   paper: 
%   
%   Global and Local Uncertainty Principles for Signals on Graphs
%
%   Authors: Nathanael Perraudin, Benjamin Ricaud, David I Shuman, Pierre
%   Vandergheynst 
%
%   ArXiv: http://arxiv.org/abs/1603.03030
%
%   Abstract of the paper
%   ---------------------
%
%   Uncertainty principles such as Heisenberg's provide limits on the
%   time-frequency concentration of a signal, and constitute an important
%   theoretical tool for designing and evaluating linear signal transforms.
%   Generalizations of such principles to the graph setting can inform
%   dictionary design for graph signals, lead to algorithms for
%   reconstructing missing information from graph signals via sparse
%   representations, and yield new graph analysis tools. While previous
%   work has focused on generalizing notions of spreads of a graph signal
%   in the vertex and graph spectral domains, our approach is to generalize
%   the methods of Lieb in order to develop uncertainty principles that
%   provide limits on the concentration of the analysis coefficients of any
%   graph signal under a dictionary transform whose atoms are jointly
%   localized in the vertex and graph spectral domains. One challenge we
%   highlight is that due to the inhomogeneity of the underlying graph data
%   domain, the local structure in a single small region of the graph can
%   drastically affect the uncertainty bounds for signals concentrated in
%   different regions of the graph, limiting the information provided by
%   global uncertainty principles. Accordingly, we suggest a new way to
%   incorporate a notion of locality, and develop local uncertainty
%   principles that bound the concentration of the analysis coefficients of
%   each atom of a localized graph spectral filter frame in terms of
%   quantities that depend on the local structure of the graph around the
%   center vertex of the given atom. Finally, we demonstrate how our
%   proposed local uncertainty measures can improve the random sampling of
%   graph signals.
%

% Author: Nathanael Perraudin
% Date  : 17 June 2014

%% Initialization
clear;
close all;

gsp_reset_seed(2);

%% Plotparameter

global SAVE ;

paramplot.save = SAVE;
paramplot.position = [100 100 300 200];
paramplot.savefig = 0;
paramplot.titlesize = 20;

%% Parameters
N = 64; % Size of the graph
Nf = N/4; % Number of filter
overlap = floor(sqrt(Nf)/2)*2; % Overlap of the filter


%% GRAPHS

G = cell(8,1);

G{1}= gsp_ring(N);
G{2} = gsp_sensor(N);
G{3} = gsp_random_regular(N);
G{4} = gsp_erdos_renyi(N,1/N*6);
G{5} = gsp_comet(N,floor(N/2));
G{6} = gsp_path(N);
W = ones(N-1,1);
W(1) = 0.1;
G{7} = gsp_modified_path(W);
W(1) = 0.01;
G{8} = gsp_modified_path(W);


% Number of graphs
Ng = length(G);

%% Compute the Fourier basis

for ii = 1:Ng 
    G{ii} = gsp_compute_fourier_basis(G{ii});
end

%% Create filters

g1 = cell(Ng,1);
g2 = cell(Ng,1);
g3 = cell(Ng,1);
g4 = cell(Ng,1);

param_filter.overlap = overlap;
for ii = 1:Ng 
    g1{ii} = gsp_design_itersine(G{ii},Nf,param_filter);
    param_filter.log = 1;
    param_filter.filter = 'itersine';
    param_filter.warping_type = 'custom';
    param_filter.warp_function = @(x) x;
    g3{ii} = gsp_design_warped_translates(G{ii},Nf,param_filter);

end

tol=1e-8;
for ii = 1:Ng 
    [unique_E,unique_E_inds]=unique(round(G{ii}.e*1/tol)*tol);
    param_filter2.approx_spectrum.x=unique_E;
    param_filter2.approx_spectrum.y=(unique_E_inds-1)/(G{ii}.N-1);  
    param_filter2.interpolation_type='monocubic';
    param_filter2.warping_type = 'spectrum_interpolation';
    param_filter2.filter = 'itersine';
    param_filter2.overlap = overlap;
    param_filter2.log = 0;
    g2{ii}  = gsp_design_warped_translates(G{ii}, Nf,param_filter2); 
    param_filter2.log = 1;
    g4{ii}  = gsp_design_warped_translates(G{ii}, Nf,param_filter2); 
end

%% Plot 2 filters
figure
gsp_plot_filter(G{2},g1{2});
title('(a) Gabor');
gsp_plotfig('FB_translates_uniform',paramplot);

figure
gsp_plot_filter(G{2},g2{2});
title('(b) Spectrum adapted Gabor');
gsp_plotfig('FB_translates_adapted',paramplot);

figure
gsp_plot_filter(G{2},g3{2});
title('(c) Wavelets');
gsp_plotfig('FB_wavelets_uniform',paramplot);

figure
gsp_plot_filter(G{2},g4{2});
title('(d) Spectrum adapted wavelets');
gsp_plotfig('FB_wavelets_adapted',paramplot);

%% Compute the uncertainty

bound1 = zeros(Ng,1);
bound2 = zeros(Ng,1);
bound12 = zeros(Ng,1);
bound22 = zeros(Ng,1);

bound3 = zeros(Ng,1);
bound4 = zeros(Ng,1);
bound32 = zeros(Ng,1);
bound42 = zeros(Ng,1);
mu = zeros(Ng,1);
for ii =1: Ng
    F = gsp_filterbank_matrix(G{ii},g1{ii});
    bound1(ii) = max(abs(F(:)))^(-1);
    bound12(ii) = 1/max(sqrt(sum(abs(F).^2,1)));
    
    F = gsp_filterbank_matrix(G{ii},g2{ii});
    bound2(ii) = max(abs(F(:)))^(-1);
    bound22(ii) = 1/max(sqrt(sum(abs(F).^2,1)));

    F = gsp_filterbank_matrix(G{ii},g3{ii});
    bound3(ii) = max(abs(F(:)))^(-1);
    bound32(ii) = 1/max(sqrt(sum(abs(F).^2,1)));

    F = gsp_filterbank_matrix(G{ii},g4{ii});
    bound4(ii) = max(abs(F(:)))^(-1);
    bound42(ii) = 1/max(sqrt(sum(abs(F).^2,1)));
    
    mu(ii) = G{ii}.mu;
end

%% Concatenate and display the result
resmat = [mu,bound12,bound22,bound32,bound42];
mat2tex(resmat,'table','%5.2f');


resmat = [mu,bound12.^(-1),bound22.^(-1),bound32.^(-1),bound42.^(-1)];
mat2tex(resmat,'table','%5.2f');
