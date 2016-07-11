%SPREAD_GABOR Display the spread of different Gabor transforms
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
%   In Figure 1, we analyze a series of signals on a random sensor network
%   of 100 vertices. Each signal is created by localizing a heat kernel to
%   be centered at vertex 1 (circled in black). To generate the four
%   different signals, we vary the value of the parameter $\tau$ in the
%   heat kernel. We plot the four localized kernels in the graph spectral
%   and vertex domains in the first two columns, respectively. The more we
%   "compress" $\hat{h}$ in the graph spectral domain (i.e. we reduce its
%   spectral spreading by increasing $\tau$), the less concentrated the
%   localized atom becomes in the vertex domain. The joint vertex-frequency
%   representation $|A_g T_1h_\tau(i,k)|$ of each signal is shown in the
%   third column, which illustrates the trade-off between concentration in
%   the vertex and the spectral domains. The concentration of these graph
%   Gabor transform coefficients is the quantity bounded by the uncertainty
%   principle presented in our uncertainty principle Theorem. In the last
%   row of the Figure 1, $\tau=\infty$ which leads to a Kronecker delta for
%   the kernel and a constant on the vertex domain. On the contrary, when
%   the kernel is constant, with $\tau=0$ (top row), the energy of the
%   graph Gabor coefficients stays concentrated around one vertex but
%   spreads along all frequencies.
%
%   .. figure::
%
%      Graph Gabor transform of four different signals
%
%      First column: the kernel $h_{\tau}(\lambda)$ is shown in red and the
%      localized kernel $\widehat{f_{\tau}}$ is shown in blue, both in the
%      graph spectral domain. Second column: the signal $f_{\tau}$ in the
%      vertex domain (the center vertex 1 is circled). Third column:
%      $|A_g T_1h_\tau (i,k)|$, the absolute value of the Gabor
%      transform coefficients for each vertex $i$ and each of the 20
%      frequency bands $k$. Fourth column: since it is hard to see where on
%      the graph the transform coefficients are concentrated when the nodes
%      are placed on a line in the third column, we display the value
%      $\sum_{k=0}^{19}  |A_g T_1h_\tau (i,k)|$ on each vertex $i$ in
%      the network. This figure illustrates the tradeoff between the vertex
%      and the frequency concentration.
%
%   References: perraudin2016global
%


% Author: Nathanael Perraudin
% Date  : 09.12.2013

%% Initialization


gsp_reset_seed(0);
clear;
close all;

%% Plotparameter

global SAVE ;
    
paramplot.save = SAVE;
paramplot.position = [100 100 300 200];
paramplot.savefig = 0;



%% Sensor graph

% Create the graph
N = 100;
Nf = round(sqrt(N))*2;
paramsensor.distribute = 1;
G = gsp_sensor(N,paramsensor);
G = gsp_compute_fourier_basis(G);
G = gsp_spectrum_cdf_approx(G);
param.filter = 'itersine';
param.overlap = 2;
param.log = 0;
g = gsp_design_itersine(G,Nf,param);
%g = gsp_design_warped_translates(G,Nf,param);
% Sorting
[~,inds] = sort(G.U(:,2),'descend');

% position
[~,tmp] = max(sum(abs(G.U')));
node = tmp;
 nodes = 1;
i0 = inds(node);
% parameter for the window 
tau = [0 4 10 1000];
b = 0;
figure(1);
for ii=1:numel(tau)
    k = @(x) exp(-(x/G.lmax-b).^2*tau(ii)^2);
    %fhat = fhat/norm(fhat);

    f = gsp_localize(G,k,i0);
    nf = norm(f);
    f = f/nf;
    fhat = gsp_gft(G,f);
    subplot(numel(tau),4,1+(ii-1)*4);
    gsp_plot_signal_spectral(G,fhat)
%     gsp_plotfig(['gaussian_kernel_spectral_',num2str(tau(ii))],paramplot);
    hold on;
    x = linspace(0, G.lmax,1000);
    plot(x,k(x)/nf,'r','LineWidth',3);
    %axis([0 G.lmax 0 1])
    subplot(numel(tau),4,2+(ii-1)*4);

    paramsum.vertex_highlight = i0;
    gsp_plot_signal(G,f,paramsum);
    caxis([0 1]);

%     gsp_plotfig(['gaussian_kernel_vertex',num2str(tau(ii))],paramplot);
    % Compute the coefficient

    Af = gsp_vec2mat(gsp_filter_analysis(G,g,f),Nf);
    % Display the results

    subplot(numel(tau),4,3+(ii-1)*4);
    gsp_plot_sgram(G,Af(inds,:)');
%     gsp_plotfig(['Agg_sensor_',num2str(tau(ii))],paramplot);
     caxis([0 0.2]);
    ylabel('Freq. band')
    
    subplot(numel(tau),4,4+(ii-1)*4);
    gsp_plot_signal(G,sum(abs(Af')));
    caxis([0 1.2])

end
%%
paramplot.position = [100 100 800 600];
gsp_plotfig('spread_graph_spectral_frame',paramplot)

%%
% figure(2)
% paramplot.show_sum = 0;
% gsp_plot_filter(G,g,paramplot);
% gsp_plotfig()
