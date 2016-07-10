%MODIFIED_PATH_EIGENVECTORS Display some eigenvectors for a modified path
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
%   In Figure 1, we display the eigenvector associated with the largest
%   graph Laplacian eigenvalue for a modified path graph of 100 nodes, for
%   several values of the weight $W_{12}$. Observe that the shape of the
%   eigenvector has a sharp local change at node 1.
%
%   .. figure::
%
%      Different graphs eigenvectors
%
%      Eigenvectors associated with the largest graph Laplacian eigenvalue
%      of the modified path graph with $100$ nodes, for different values of
%      $W_{12}$. As the distance between the first two nodes increases, the
%      eigenvector becomes sharply peaked.
%
%   References: perraudin2016global
%
  

%% Initialization
clear

close all;

%% Plotparameter
global SAVE ;

paramplot.save = SAVE;
paramplot.position = [100 100 300 200];
paramplot.savefig = 0;

%%

% Number of nodes
N=100;

% values to be tested
val_s=[1,10,100,1000];


W = ones(N-1,1);


figure
for ii=1:length(val_s);
    W(1) = 1/val_s(ii);
    
    % Create the graph
    G = gsp_modified_path(W);

    G = gsp_compute_fourier_basis(G); 

    % Save the result on the ambiguity fonction
    [~,ind] = max(max(abs(G.U)));
    subplot(round(sqrt(length(val_s))),round(sqrt(length(val_s))),ii);
    stem(G.U(:,ind));
    h = title(['$W_{12}$ = ',num2str(W(1))],'FontSize',14);
    set(h,'interpreter','latex','FontSize',14);
end
paramplot.position = [100 100 600 300];
gsp_plotfig('Aggpath_eig',paramplot);
