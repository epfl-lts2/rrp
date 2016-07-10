%GLOBAL_ILLUSTRATION Laplacian eigenvector localization
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
%   Let us consider the two manifolds (surfaces) embedded in $\Rbb^3$ and
%   shown in Figures 1 and 2. The first one is a flat square. The second is
%   identical except for the center where it contains a spike. We sample
%   both of these manifolds uniformly across the $x$-$y$ plane and create a
%   graph by connecting the $8$ nearest neighbors with weights depending on
%   the distance ($W_{ij}=e^{-d_{ij}/\sigma}$). The energy of each
%   Laplacian eigenvector of the graph arising from the first manifold is
%   not concentrated on any particular vertex; i.e.,
%   $\max_{i,\ell}|u_\ell(i)| << 1$, where $u_\ell$ is the eigenvector
%   associated with eigenvalue $\lambda_\ell$. However, the graph arising
%   from the second manifold does have a few eigenvectors, such as
%   eigenvector 3 shown in Figure 3, whose energy is highly concentrated on
%   the region of the spike; i.e: $\max_{i,\ell}|u_\ell(i)| \approx 1$.
%   Yet, the Laplacian eigenvectors of this second graph whose energy
%   resides primarily on the flatter regions of the manifold, such as
%   eigenvector 17 shown in Figure 5, are not too concentrated on any
%   single vertex. Rather, they more closely resemble some of the Laplacian
%   eigenvectors of the graph arising from the first manifold.
%
%   We discretize two different manifolds by sampling uniformly across the
%   $x$-$y$ plane. Due its bumpy central part, the graph arising from
%   manifold 2 has a graph Laplacian eigenvector (shown in the middle row
%   of the right column) that is highly concentrated in both the vertex and
%   graph spectral domains. However, the eigenvectors of this graph whose
%   energy primarily resides in the flatter parts of the manifold (such as
%   the one shown in the bottom row of the right column) are less
%   concentrated, and some closely resemble the Laplacian eigenvectors of
%   the graph arising from the flat manifold 1 (such as the corresponding
%   eigenvector shown in the bottom row of the left column.
%
%   .. figure::
%
%      Manifold 2
%
%      
%
%   .. figure::
%
%      Manifold 1
%
%      
%
%   .. figure::
%
%      Eigenvector no 3, manifold 2
%
%      
%
%   .. figure::
%
%      Eigenvector no 5, manifold 1
%
%      
%
%   .. figure::
%
%      Eigenvector no 17, manifold 2
%
%      
%
%   .. figure::
%
%      Eigenvector no 13, manifold 1
%
%      
%   
%   References: perraudin2016global
%

%% Initialization
clear
close all;

%% Plotparameter
global SAVE ;
if numel(SAVE)==0
    SAVE = 0;
end

paramplot.save = SAVE;
paramplot.savefig = 0;


%% Parameter of sampling

tau = 0.1;
N = 30;
s = 3;
    [X,Y] = meshgrid(-s:2*s/N:s,-s:2*s/N:s);
Z = max(1./(tau+ X.^2 + Y.^2),1);

Nf = 200;
[Xf,Yf] = meshgrid(-s:2*s/Nf:s,-s:2*s/Nf:s);
Zf = max(1./(tau+ Xf.^2 + Yf.^2),1);

%% Create the graphs

P = [X(:),Y(:),Z(:)];
paramnn.k = 8;
paramnn.sigma = 0.03;
paramnn.rescale = 1;
G = gsp_nn_graph(P,paramnn);
Gbase = gsp_nn_graph([X(:),Y(:)],paramnn);

% Compute Fourier basis
G = gsp_compute_fourier_basis(G);
Gbase = gsp_compute_fourier_basis(Gbase);



%% Plot the results

[val,ind] = max(abs(G.U));
[~,ind2] = max(val);
ind1 = ind(ind2);

paramplot.position = [100 100 400 300];

figure(1)
hold off
surface(Xf,Yf,Zf,-Zf,'EdgeColor','none')
view([-25,39])
caxis([min(-Zf(:)),max(-Zf(:))]);
hold on
scatter3(X(:),Y(:),Z(:),10*ones(size(X(:))),'r','filled');
axis([min(X(:)) max(X(:)) min(Y(:)) max(Y(:)) min(Z(:)) max(Z(:))])
if ~SAVE
    title('Manifold 2')
end
gsp_plotfig('original_manifold_1',paramplot)

figure(2)
hold off
surface(Xf,Yf,ones(size(Xf)),-ones(size(Xf)),'EdgeColor','none')
view([-25,39])
caxis([min(-Zf(:)),max(-Zf(:))]);
hold on
scatter3(X(:),Y(:),ones(size(X(:))),10*ones(size(X(:))),'r','filled');
axis([min(X(:)) max(X(:)) min(Y(:)) max(Y(:)) min(Z(:)) max(Z(:))+eps])
if ~SAVE
    title('Manifold 1')
end
gsp_plotfig('original_manifold_2',paramplot)

paramplot.position = [100 100 300 200];

figure(3)
imagesc(reshape(G.U(:,ind2),size(X)))
if ~SAVE
    title(['Eigenvector no ', num2str(ind2)]);
end
colorbar
axis off
gsp_plotfig('manifold_1_eig_1',paramplot)


figure(4)
imagesc(reshape(-Gbase.U(:,5),size(X)))
if ~SAVE
    title('Eigenvector no 5');
end
colorbar
axis off

gsp_plotfig('manifold_2_eig_1',paramplot)

figure(5)
imagesc(reshape(G.U(:,17),size(X)))
if ~SAVE
    title('Eigenvector no 17');
end
colorbar
axis off
gsp_plotfig('manifold_1_eig_2',paramplot)

figure(6)
imagesc(reshape(-Gbase.U(:,13),size(X)))
if ~SAVE
    title('Eigenvector no 13');
end
colorbar
axis off

gsp_plotfig('manifold_2_eig_2',paramplot)