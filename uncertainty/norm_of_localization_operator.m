%NORM_OF_LOCALIZATION_OPERATOR This script illustrates a property of the localization operator
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
%   This experiment
%   ---------------
%
%   This experiment illustrates the effect of the graph
%   structure on the norms of localized functions.
% 
%   We take the kernel to be localized to be a heat kernel of the form
%   $\hat{g}(\lambda_\ell) = e^{-\tau \lambda_\ell}$, for some constant
%   $\tau>0$. We localize the kernel $\hat{g}$ to be centered at each
%   vertex $i$ of the graph with the operator $T_i$, and we compute and
%   plot their $\ell^2$-norms $\|T_ig\|_2$. The figure shows that when a
%   center node $i$ and its surrounding vertices are relatively weakly
%   connected, the $\ell^2$-norm of the localized heat kernel is large, and
%   when the nodes are relatively well connected, the norm is smaller.
%   Therefore, the norm of the localized heat kernel may be seen as a
%   measure of vertex centrality.
%
%   .. figure::
%
%      Random sensor
%
%      
%
%   .. figure::
%
%      Comet
%
%      
%
%   .. figure::
%
%      Grid
%
%      
%
%   .. figure::
%
%      Airfoil
%
%      
%
%   .. figure::
%
%      Minnesota
%
%
%
%   .. figure::
%
%      The heat kernel
%
%      $\hat{g}(\lambda_\ell)=e^{-10\frac{\lambda_\ell}{\lambda_\text{max}}}$
%   
%   References: perraudin2016global
%

%% Initialisation

clear;
close all;



%% Parameters
global SAVE
paramplot.show_edges = 1;
paramplot.save = SAVE;
paramplot.savefig = 0;
paramplot.titleweight = 'bold';
paramplot.position = [100,100, 450, 300];
% Kernel
tau = 1;

g = @(x) exp(-tau*x);

%% Create graphs
G1 = gsp_sensor();
G2 = gsp_comet(16,8);
G3 = gsp_2dgrid(16);
G3.vertex_size = 50;
G4 = gsp_airfoil();
G5 = gsp_minnesota();
G6 = gsp_swiss_roll();




%% Translate the kernel
G1 = gsp_estimate_lmax(G1);
s1 = sqrt(G1.N)*gsp_filter_analysis(G1,g,eye(G1.N));

G2 = gsp_estimate_lmax(G2);
s2 = sqrt(G2.N)*gsp_filter_analysis(G2,g,eye(G2.N));

G3 = gsp_estimate_lmax(G3);
s3 = sqrt(G3.N)*gsp_filter_analysis(G3,g,eye(G3.N));

G4 = gsp_estimate_lmax(G4);
s4 = sqrt(G4.N)*gsp_filter_analysis(G4,g,eye(G4.N));

G5 = gsp_estimate_lmax(G5);
s5 = sqrt(G5.N)*gsp_filter_analysis(G5,g,eye(G5.N));


G6 = gsp_estimate_lmax(G6);
s6 = sqrt(G6.N)*gsp_filter_analysis(G6,g,eye(G6.N));

%% Compute the norm and plot the results

% n1 = (sum(abs(s1),1)./sqrt(sum(abs(s1).^2,1)))';
% n2 = (sum(abs(s2),1)./sqrt(sum(abs(s2).^2,1)))';
% n3 = (sum(abs(s3),1)./sqrt(sum(abs(s3).^2,1)))';
% n4 = (sum(abs(s4),1)./sqrt(sum(abs(s4).^2,1)))';
% n5 = (sum(abs(s5),1)./sqrt(sum(abs(s5).^2,1)))';

n1 = sqrt(sum(abs(s1).^2,1))';
n2 = sqrt(sum(abs(s2).^2,1))';
n3 = sqrt(sum(abs(s3).^2,1))';
n4 = sqrt(sum(abs(s4).^2,1))';
n5 = sqrt(sum(abs(s5).^2,1))';
n6 = sqrt(sum(abs(s6).^2,1))';


%%
figure()
gsp_plot_signal(G1,n1,paramplot);
title('Random sensor');
gsp_plotfig('heat_translate_sensor',paramplot);

figure()
gsp_plot_signal(G2,n2,paramplot);
title('Comet');

gsp_plotfig('heat_translate_comet',paramplot);

figure()
gsp_plot_signal(G3,n3,paramplot);
title('Grid');
gsp_plotfig('heat_translate_grid',paramplot);

figure()
gsp_plot_signal(G4,n4,paramplot);
title('Airfoil')
gsp_plotfig('heat_translate_airfoil',paramplot);
%%
figure()
gsp_plot_signal(G5,n5,paramplot);
title('Minnesota')
gsp_plotfig('heat_translate_minnesota',paramplot);


%%
figure()
paramplot.show_sum = 0;
gsp_plot_filter(G3,g,paramplot)
title('Heat kernel')
gsp_plotfig('heat_translate_kernel',paramplot);

