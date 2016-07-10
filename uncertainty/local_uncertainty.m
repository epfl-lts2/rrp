%LOCAL_UNCERTAINTY Local unceratainty illustration
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
%   Let us concentrate on the case where $p=\infty$ of the Local
%   uncertainty Theorem. It tells us that the concentration of
%   $A_{g}T_{i_0}g_{k_0}$ is limited by
%   $\frac{1}{\|T_{\tilde{i}}g_{\tilde{k}_{i_0,k_0}}\|_{2}}$. One question
%   is to what extent this quantity is local or reflects the local behavior
%   of the graph. As a general illustration for this discussion, we present
%   in the figures bellow quantities related to the local uncertainty of a
%   random sensor network of $100$ nodes evaluated for two different values
%   of $k$ (one in each column) and all nodes $i$.
%
%   .. figure::
%
%      Gabor filterbank used
%
%      
%
%   .. figure::
%
%      Local sparsity level
%
%      
%
%   .. figure::
%
%      Theorem bound
%
%      
%
%   .. figure::
%
%      Hope distanz of between $i_0$ and $\tilde{i}$
%
%      
%
%   .. figure::
%
%      Approximated bound
%
%      
%
%   .. figure::
%
%      $\tilde{k}$
%
%      
%
%   .. figure::
%
%      Relative difference between bound approximation and Theorem bound
%
%      
%   
%   References: perraudin2016global
%



%% Initialisation
clear;
close all;

global SAVE;
global sf
% plotting paramter
paramplot.save = SAVE;
paramplot.position = [100 100 300 200];
paramplot.savefig = 0;

gsp_start;
gsp_reset_seed(0);

%% Prepare the graph

N = 100;
Nf = 7;
if ~numel(sf)
    sf = 0;
end

paramgraph.distribute = 0;

G = gsp_sensor_old(N,paramgraph);
% G = gsp_2dgrid(10);
G = gsp_compute_fourier_basis(G);

%% Design the gabor filterbank

g = gsp_design_itersine(G,Nf);
% Frame bounds
LB = gsp_filterbank_bounds(G,g);

%%

sp1 = zeros(G.N,1);
approx_upper_bound = zeros(G.N,1);
B5 = zeros(G.N,1);
theo_bound = zeros(G.N,1);
ktilde = zeros(G.N,1);
itilde = zeros(G.N,1);
hop_dist = zeros(G.N,1);

for ii = 1:G.N
    
    gi = gsp_localize(G,g{sf+1},ii);
    g2i = gsp_localize(G,@(x) (g{sf+1}(x)).^2,ii);
    Aggi = gsp_filter_analysis(G,g,gi);
    
    [max_Aggi,k] = max(abs(gsp_vec2mat(Aggi,Nf)'));
    [max_Aggi, itilde(ii)] = max(max_Aggi);
    ktilde(ii) = k(itilde(ii));
    hop_dist(ii) = gsp_hop_distanz(G,ii,itilde(ii));
    sp1(ii) = (norm(Aggi(:),2)/max_Aggi)^(-1);
    approx_upper_bound(ii) = (sqrt(LB*G.N)/norm(gi,2))^(-1); 
    
    gki = gsp_localize(G,g{ktilde(ii)},itilde(ii));
    theo_bound(ii) = (sqrt(LB*G.N)/norm(gki,2))^(-1); 


end

B3 = 1./max(abs(G.U),[],2);


%%
figure
gsp_plot_filter(G,g);
%title('Gabor filterbank')
gsp_plotfig('LU_filters',paramplot)

%%
figure()

gsp_plot_signal(G,theo_bound)
%title('Theorem bound');
caxis([0, max(sp1)]);
gsp_plotfig(['LU_theo_bound_1',num2str(sf+1)],paramplot);


%%
figure;
gsp_plot_signal(G,sp1);
%h= title('Local uncertainty for filter $\frac{\|B_{g}T_{i0}g_{k0}\|_{2}}{||B_{g}T_{i0}g_{k0}\|_{\infty}}$');
%set(h,'interpreter','latex','FontSize',14);
%title(['Local uncertainty for filter ',num2str(sf)]);
caxis([0, max(sp1)]);
gsp_plotfig(['LU_local_uncertainty',num2str(sf+1)],paramplot);



%%

figure

gsp_plot_signal(G,hop_dist)
%h = title('Hope distanz of between $i_0$ and $\tilde{i}$');
%set(h,'interpreter','latex','FontSize',14);
gsp_plotfig(['LU_hope_dist_1',num2str(sf+1)],paramplot);

%%
figure
gsp_plot_signal(G,approx_upper_bound);
% h= title('Approximated bound $\frac{\sqrt{A}}{\|T_{i0}g_{k0}\|_2}$');
% set(h,'interpreter','latex','FontSize',14);
%title(['Approximated bound for filter ',num2str(sf)]);
caxis([0, max(sp1)]);
gsp_plotfig(['LU_localised_g0',num2str(sf+1)],paramplot);

%%
% figure
% sig = abs(B1-B5)./B1;
% sig(abs(sig)<(1e-14)) = 0;
% gsp_plot_signal(G,sig);
% h = title(['Rel. err. between $\| T_i g_{k_0}^2 \|_\infty$ and $\max_k \| T_i g_{k_0} g_k \|_\infty$, filter ',num2str(sf)]);
% set(h,'interpreter','latex','FontSize',14);
% gsp_plotfig(['LU_rel_error_1',num2str(sf)],paramplot);
figure
sig = ktilde-1;
gsp_plot_signal(G,sig);
caxis([0 2])
%h = title(['$\tilde{k}(i)$, filter ',num2str(sf)]);
%set(h,'interpreter','latex','FontSize',14);
gsp_plotfig(['LU_rel_error_1',num2str(sf+1)],paramplot);


%%

figure
gsp_plot_signal(G,abs(sp1-approx_upper_bound)./sp1);
%h = title(['Rel. err. between $\| T_i g_{k_0} \|_2^2$ and $\max_k \| T_i g_{k_0} g_k \|_\infty$, filter ',num2str(sf)]);
%set(h,'interpreter','latex','FontSize',14);
gsp_plotfig(['LU_rel_error_2',num2str(sf+1)],paramplot);


