%MODIFIED_PATH_GABOR Show global and local bound on the modified path
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
%   On a $64$ node modified path graph for details), we compute the graph
%   Gabor transform of the signals $f_1=T_1g_0$ and $f_2=T_{64}g_0$. In
%   Figure \ref{fig:path_node_away_ambiguity}, we show the evolution of the
%   graph Gabor transforms of the two signals with respect to the distance
%   $d=1/W_{12}$ from the first to the second vertex in the graph. As the
%   first node is pulled away, a localized eigenvector appears centered on
%   the isolated vertex. Because of this, as this distance increases, the
%   signal $f_1$ becomes concentrated in both the vertex and graph spectral
%   domains, leading to graph Gabor transform coefficients that are highly
%   concentrated. However, since the graph modification is local, it does
%   not drastically affect the graph Gabor transform coefficients of the
%   signal $f_2$, whose energy is concentrated on the far end of the path
%   graph.
%
%   In Figure 1, we plot the evolution of the uncertainty bounds as well as
%   the concentration of the Gabor transform coefficients of $f_1$ and
%   $f_2$. The global uncertainty bound from our Theorem tells us that
%   $s_1(A_g f)\leq\max_{i,k}||T_i g_k||_2,\hbox{ for any signal }f.$
%   The local uncertainty bound from our Theorem tells us that 
%   $s_1(\A_\gg T_{i_0}g_{k_0})\leq||T_{\tilde{i}_{i_0,k_0}}g_{\tilde{k}_{i_0,k_0}}||_2,\hbox{ for all }i_0\hbox{ and }k_0.$
%   Thus, we can view the global uncertainty bound as an upper bound on all 
%   of the local uncertainty bounds. In fact the bumps in the global 
%   uncertainty bound in Figure~\ref{fig:bound Agf} correspond to the local 
%   bound with $i_0=1$ and   different frequency bands $k_0$. 
%   We plot the local bounds for $i_0=1$ and $k_0=0$ and $k_0=2$. 
%
%   .. figure::
%
%      Evolution of the bounds
%
%      Concentration of the graph Gabor coefficients of $f_1=T_1g_0$ and
%      $f_2=T_{64}g_0$ with respect to the distance between the first two
%      vertices in the modified path graph, along with the upper bounds on
%      this concentration from  Theorem  \ref{Co:Lieblocgraph} (global
%      uncertainty) and Theorem \ref{theo:local_uncertainty} (local
%      uncertainty). Each bump of the global bound corresponds to a local
%      bound of a given spectral band of node $1$. For clarity, we plot
%      only bands $\widehat{g_0}$ and $\widehat{g_2}$ for node $1$. For
%      node $64$, the local bound is barely affected by the change in graph
%      structure, and the sparsity levels of the graph Gabor transform
%      coefficients of $T_{64}g_0$ also do not change much.
%   
%   References: perraudin2016global
%

% Author: Nathanael Perraudin
% Date  : 09.12.2015

%% Initialization
clear

close all;

%% Plotparameter
global SAVE ;

paramplot.save = SAVE;
paramplot.position = [100 100 300 200];
paramplot.savefig = 0;

%% Parameter
% Number of node of the graph
N=64;

% Norm chosen
p = 1;

k = 0; %filter number

% values to be tested
val=9.^(0:-0.05:-4);
dist1node = 1./val;

%% graph
W = ones(N-1,1);
G = gsp_modified_path(W);
G = gsp_compute_fourier_basis(G);
G = gsp_spectrum_cdf_approx(G);

%% Create the filter-bank
% Warping to recover the gabor case
param.filter = 'itersine';
param.log = 0;
Nf = 16;
g = gsp_design_warped_translates(G,Nf,param);
%g = gsp_design_itersine(G,Nf,param);
% Drop uncomplete filters

%% Loop
sp1 = zeros(size(val));
sp7 = zeros(size(val));
sp10 = zeros(size(val));
gb = zeros(size(val));
lb1 = zeros(size(val));
lb13 = zeros(size(val));
lb7 = zeros(size(val));
lb10 = zeros(size(val));
mu = zeros(size(val));

% printii =(mod(1:length(val),10)==1) .*((1:length(val))<50 );
printii = zeros(1,length(val));
printii([1,23, 27 31, 41]) = 1;
for ii=1:length(val);


    % Create the graph
    W = ones(N-1,1);
    W(1) = val(ii);
    G = gsp_modified_path(W);
    G = gsp_compute_fourier_basis(G);

    % Signals
    f1=zeros(N,1);
    f1(1)=1;
    f1 = gsp_filter_analysis(G,g{k+1},f1);

    f7=zeros(N,1);
    f7(round(N/2))=1;    
    f7 = gsp_filter_analysis(G,g{k+1},f7);
    
    f10 = zeros(N,1);
    f10(N) = 1;
    f10 = gsp_filter_analysis(G,g{k+1},f10);

    mu(ii) = G.mu;   

    
    
    Af1 = gsp_vec2mat(gsp_filter_analysis(G,g,f1),Nf);
    Af7 = gsp_vec2mat(gsp_filter_analysis(G,g,f7),Nf);
    Af10 = gsp_vec2mat(gsp_filter_analysis(G,g,f10),Nf);
    
    sp1(ii) = norm(Af1(:),2)/norm(Af1(:),p);
    sp7(ii) = norm(Af7(:),2)/norm(Af7(:),p);
    sp10(ii) = norm(Af10(:),2)/norm(Af10(:),p);
    
    % Upperbound  
    gb(ii) = sp_spectral_frame_bound(G,g,p);
    lb1(ii) = sp_spectral_frame_local_bound(G,g,p,1,k+1);
    lb13(ii) = sp_spectral_frame_local_bound(G,g,p,1,3);
    lb7(ii) = sp_spectral_frame_local_bound(G,g,p,round(N/2),k+1);
    lb10(ii) = max(sp_spectral_frame_local_bound(G,g,p,N,1:length(g)));
    

    % Plot the gabor transform
    if printii(ii)
       figure()
       gsp_plot_sgram(G,abs(Af1') );
    %   title(['d = ',num2str(1/val(ii))]);
%         caxis([0,max(abs(Af1(:)))])

       ylabel('Freq. band');
       gsp_plotfig(['2Aggpath_',num2str(1/val(ii))],paramplot);


       figure()
       gsp_plot_sgram(G,abs(Af10'));
       ylabel('Freq. band');
%         caxis([0,max(abs(Af1(:)))])
     %  title(['node 10, d = ',num2str(1/val(ii))]);
       gsp_plotfig(['2Aggpath10_',num2str(1/val(ii))],paramplot);
       
       figure()
       gsp_plot_signal(G,f1);
%         caxis([0,max(abs(Af1(:)))])
     %  title(['node 10, d = ',num2str(1/val(ii))]);
       gsp_plotfig(['2Aggpath_sig_',num2str(1/val(ii))],paramplot);
       
    end
    

end




%% Plot the evolution of the bounds
figure(1)
semilogx(dist1node,sp1,'b',...
    dist1node,sp10,'r',...
    dist1node,gb,'-kx',...
    dist1node,lb1,'-bo',...
    dist1node,lb13,'-go', ...
    dist1node,lb10,'-ro');
% h = legend( ' $\frac{\| \mathcal{A}_g \delta_1 \|_2}{\| \mathcal{A}_g \delta_{1} \|_1}$ $~$', ...
%      '$\frac{\| \mathcal{A}_g \delta_{64} \|_2}{\| \mathcal{A}_g \delta_{64} \|_1}$ $~$', ...
%     'Global bound $~$ $~$',...
%     'Local bound: node 1, band $g_0$ $~$',...
%     'Local bound: node 1, band $g_2$ $~$',...
%     'Local bound: node 64, all bands$~$ $~$',...
%     'Location','East');
h = legend( ' $\frac{\| \mathcal{A}_g T_1g_0  \|_2}{\| \mathcal{A}_g T_1g_0 \|_1}$ $~$', ...
     '$\frac{\| \mathcal{A}_g T_{64}g_0 \|_2}{\| \mathcal{A}_g T_{64}g_0  \|_1}$ $~$', ...
    'Global bound $~$ $~$',...
    'Local bound: node 1, band $\hat{g}_0$ $~$',...
    'Local bound: node 1, band $\hat{g}_2$ $~$',...
    'Local bound: node 64, all bands$~$ $~$',...
    'Location','East');
set(h,'interpreter','latex','FontSize',14,'Position', [0.5981 0.4425 0.4237 0.2951]);
xlabel('Distance between nodes 1 and 2');
set( findobj(gca,'type','line'), 'LineWidth', 2);
axis([min(dist1node) max(dist1node) 0 1])
gsp_plotfig('2BoundAgg21_path' );



