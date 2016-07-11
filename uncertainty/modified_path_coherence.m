%MODIFIED_PATH_COHERENCE Illustrates how the relation between coherence and uncertainty
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
%   Experiment 1
%   ------------
%
%   We start with a standard path graph of $10$ nodes equally spaced (all
%   edge weights are equal to one) and we move the first node out to the
%   left; i.e., we reduce the weight between the first two nodes (see
%   Figure 1). The weight is related to the distance by $W_{12}=1/d(1,2)$
%   with $d(1,2)$ being the distance between nodes 1 and 2. When the weight
%   between nodes 1 and 2 decreases, the eigenvector associated with the
%   largest eigenvalue of the Laplacian becomes more concentrated, which
%   increases the coherence $\mu_{G}$ (Figure 2).
%
%   .. figure::
%
%      Example of a modified path graph with $10$ nodes
%
%      
%
%   .. figure::
%
%      Graph coherence
%
%      Evolution of the coherence of the modified path graph with respect
%      to the distance between nodes $1$ and $2$. As the degree of the
%      comet's center vertex increases or the first node of the modified
%      path is pulled away, the coherence $\mu_{G}$ tends to the limit
%      value $\sqrt{\frac{N-1}{N}}$
%
%   Experiment 2
%   ------------
%
%   .. .||f||_p ||\hat{f}||_p \geq \mu_{G}^{1-\frac{2}{p}}||f||_{2}^2  
%   ..  for p in [1 2] (1)
%
%   .. math:: \|f\|_{p}\|\hat{f}\|_{p}\geq \mu_{G}^{1-\frac{2}{p}}\|f\|_{2}^2 ,\qquad p\in[1,2], \hspace{1cm} (1)
%
%   Figure 3 shows the computation of the quantities involved in (1), with
%   $p=1$ and different graph's taken to be the modified path graphs of the
%   modified path, with different distances between the first two vertices.
%   We show the lefthand side of (1) for two different Kronecker deltas,
%   one centered at vertex 1, and one centered at vertex 10. We have seen
%   in |global_illustration| that as the distance between the first two
%   vertices increases, the coherence increases, and therefore the lower
%   bound on the right-hand side of (1) decreases. For $\delta_1$, the
%   uncertainty quantity on the left-hand side of (1) follows a similar
%   pattern. The intuition behind this is that as the weight between the
%   first two vertices decreases, a few of the eigenvectors start to have
%   local jumps around the first vertex (see |modified_path_eigenvectors|
%   ). As a result, we can sparsely represent $\delta_1$ as a linear
%   combination of those eigenvectors and $||\widehat{\delta_1}||_1$ is
%   reduced. However,  since there are not any eigenvectors that are
%   localized around the last vertex in the path graph, we cannot find a
%   sparse linear combination of the graph Laplacian eigenvectors to
%   represent $\delta_{10}$. Therefore, its uncertainty quantity on the
%   left-hand side of (1) does not follow the behavior of the lower bound.
%
%   .. figure::
%
%      Evolution of bounds
%
%      Numerical illustration of the $l^p$-norm uncertainty principle on a
%      sequence of modified path graphs with different mutual coherences
%      between the canonical basis of deltas and the graph Laplacian
%      eigenvectors. For each modified path graph, the weight $W_{12}$ of
%      the edge between the first two vertices is the reciprocal of the
%      distance shown on the horizontal axis. The black crosses show the
%      lower bound on the right-hand side of (1), with $p=1$. The blue and
%      red lines show the corresponding uncertainty quantity on the
%      left-hand side of  \eqref{eq:rictorr}, for the graph signals
%      $\delta_1$ and $\delta_{10}$, respectively.
%
%   Experiment 3
%   ------------
%
%   The Hausdorff-Young Theorem for graph signals is stated as:
%
%   Let $\mu_{G}$ be the coherence between the graph Fourier and canonical
%   bases of a graph $G$. Let $p,q>0$ be such that
%   $\frac{1}{p}+\frac{1}{q}=1$. For any signal $f \in \mathbb{C}^N$
%   defined on $G$ and $1 \leq p \leq 2$, we have
%   
%   .. .|| \hat f ||_q \leq \mu_{G}^{1-\frac{2}{q}} || f ||_p, (2)
%
%   .. math:: \| \hat f \|_q \leq \mu_{G}^{1-\frac{2}{q}} \| f \|_p, \hspace{1cm} (2)
%   
%   and conversly:
%   
%   .. .|| \hat f \|_q \geq \mu_{G}^{1-\frac{2}{q}} || f ||_p. (3)
%
%   .. math:: \| \hat f \|_q \geq \mu_{G}^{1-\frac{2}{q}} \| f \|_p. \hspace{1cm} (3)
%
%
%   Continuing with the modified path graphs, we illustrate the bounds of
%   the Hausdorff-Young Hausdorff-Young Theorem. this example, we take the
%   signal $f$ to be $\delta_1$, a Kronecker delta centered on the first
%   node of the modified path graph. As a consequence, $\|\delta_1\|_p=1$
%   for all $p$, which makes it easier to compare the quantities involved
%   in the inequalities. For this example, the bounds of Hausdorff-Young
%   Theorem are fairly close to the actual values of
%   $\|\hat{\delta_1}\|_q$.
%
%   .. figure::
%
%      Illustration of the bounds of the Hausdorff-Young inequalities
%
%      Illustration of the bounds of the Hausdorff-Young inequalities for
%      graph signals on the modified path graphs with $f=\delta_1$. The
%      quantities in (2) and (3) for $q=1,\frac{4}{3},4,$ and $\infty$.
%
%   .. figure::
%
%      Illustration of the bounds of the Hausdorff-Young inequalities
%
%      The same quantities with respect of the sparsity level
%   
%   References: perraudin2016global
%



% Author: Nathanael Perraudin
% Date  : 

%% Initialization
clear

close all;

%% Plotparameter
global SAVE ;

paramplot.save = SAVE;
paramplot.position = [100 100 300 200];
paramplot.savefig = 0;

paramplot2 = paramplot;
paramplot2.position = [100 100 450 300];

paramplot3 = paramplot;
paramplot3.position = [100 100 300 200];

paramplot4 = paramplot;
paramplot4.position = [100 100 300 100];
%% graph

% Number of nodes
N=10;

% Two signals
f=zeros(N,1);
f(1)=1;

f10 = zeros(N,1);
f10(10) = 1;


% values of distances to be tested
val=9.^(0:-0.05:-3);
dist1node = 1./val;

% Prepare variable for saving the data
savedat=zeros(size(val,2),N);
n1dat=zeros(size(val));
n1dat10=zeros(size(val));
n13dat=zeros(size(val));
n2dat=zeros(size(val));
n4dat=zeros(size(val));
ninfdat=zeros(size(val));
mu = zeros(size(val));


for ii=1:length(val);

    
    % Create the graph
    W = ones(N-1,1);
    W(1) = val(ii);
    G = gsp_modified_path(W);

    G = gsp_compute_fourier_basis(G);
    
    % Save the maximum of the eigenvectors
    [savedat(ii,:),~]=max(abs(G.U));
    
    % Save the norm 1 of the Fourier transform
    fhat=gsp_gft(G,f);
    fhat10=gsp_gft(G,f10);
    n1dat(ii) = norm(fhat,1);
    n1dat10(ii) = norm(fhat10,1);
    
    % Save other norms of the Fourier transform
    n13dat(ii) = norm(fhat,4/3);
    n2dat(ii) = norm(fhat,2);
    n4dat(ii) = norm(fhat,4);
    ninfdat(ii) = max(abs(fhat(:)));

    % Save the result on the ambiguity fonction
    [mu(ii),~] = max(savedat(ii,:));

    % Plot a example of the modifed path graph
    if ii==30
        %%
        figure;
        gsp_plot_signal(G,max(abs(G.U),[],2));
        hold on;
        %quiver(0.6,0.5,-0.1,0,'k','filled','LineWidth',3,'MaxHeadSize',0.6);
        annotation('arrow',[0.25,0.1],[0.65 0.65],'Color','k','LineWidth',3);
        colorbar off
        title('Modified path graph')
        gsp_plotfig('path_node_away',paramplot4);
    end

end



%% Coherence with respect of the distance between node 1 and node 2
figure
semilogx(dist1node,mu);
xlabel('Distance between nodes 1 and 2','interpreter','latex','FontSize',14);
ylabel('$\mu_{\mathcal{G}}$: coherence', 'interpreter','latex','FontSize',14)

gsp_plotfig('mu_path',paramplot3);

%% Vertex/Fourier pure concentration
figure
semilogx(dist1node,n1dat,'b',...
    dist1node,n1dat10,'r',...
    dist1node,1./mu,'kx');
h = legend(' $\| \delta_1 \|_1  \cdot  \|\hat{\delta}_1\|_1 ~ $ ',...
    ' $\| \delta_{10} \|_1 \cdot   \|\hat{\delta}_{10}\|_1 ~ $  ',...
    'Lower bound: $\frac{1}{\mu_G} ~ $  $~$','Location','Best');
set(h,'Position',[0.4083 0.5436 0.5385 0.3024]);
set(h,'interpreter','latex','FontSize',14);
xlabel('Distance between nodes 1 and 2');
set( findobj(gca,'type','line'), 'LineWidth', 2);
gsp_plotfig('Bound11mu_path',paramplot2);



%% Hausdorf Young inequality
figure
semilogx(    dist1node,n2dat,'k',...
    dist1node,n1dat,'b',...
    dist1node,1./mu,'bx',...
    dist1node,n13dat,'c',...
    dist1node,1./sqrt(mu),'cx',...
    dist1node,n4dat,'m',...
    dist1node,mu,'rx',...
    dist1node,ninfdat,'r',...
    dist1node,sqrt(mu),'mx')
axis([1 1000 0 4])
h = legend(    'Energy: $ \|\delta_1\|_2 ~ $ ',...
    '$\|\hat{\delta_1}\|_1 ~ $ ',...
    'Lower bound on $\|\hat{\delta_1}\|_1$: $\frac{1}{\mu_G} ~ $',...
    '$\|\hat{\delta_1}\|_{4/3} ~ $ ',...    
    'Lower bound on $\|\hat{f}\|_{4/3}$: $\frac{\delta_1}{\sqrt{\mu_G}} ~$ ',...
    '$\|\hat{\delta_1}\|_4 ~ $ ',...
    'Upper bound on $\|\hat{\delta_1}\|_4$: $\sqrt{\mu_G}$',...
    '$\|\hat{\delta_1}\|_\infty ~ $ ',...
    'Upper bound on $\|\hat{\delta_1}\|_\infty$: $\mu_G ~ $');
set(h,'Position',[0.4706 0.4040 0.4939 0.5960]);
set(h,'interpreter','latex','FontSize',14);
xlabel('Distance between nodes 1 and 2','interpreter','latex','FontSize',14);
set( findobj(gca,'type','line'), 'LineWidth', 2);


gsp_plotfig('Bound11_HJ',paramplot2);



%% Hausdorf Young expressed with sparsity levels.

spsq1 = n2dat./n1dat;
spsq13 = n2dat./n13dat;
spsq2 = n2dat./n2dat;
spsq4 = n4dat./n2dat;
spsqinf = ninfdat./n2dat;

figure
semilogx(dist1node,spsq1,'m',...
    dist1node,spsqinf,'r',...
    dist1node,mu,'kx',...
    dist1node,spsq13,'b',...
    dist1node,spsq4,'c',...
    dist1node,sqrt(mu),'gx');
h = legend('$s_\infty(\delta_1)  \cdot  s_1(\hat{\delta}_1) ~ $ ',...
    '$s_1(\delta_1)  \cdot  s_\infty(\hat{\delta}_1) ~ $ ',...
    'Upper bound: $\mu_G ~ $  $~$',...
    '$s_4(\delta_1)  \cdot  s_{4/3}(\hat{\delta}_1) ~ $ ',...
    '$s_{4/3}(\delta_1)  \cdot  s_4(\hat{\delta}_1) ~ $ ',...
    'Upper bound: $\sqrt{\mu_G} ~ $  $~$','Location','Best');
set(h,'interpreter','latex','FontSize',14);
xlabel('Distance between nodes 1 and 2','interpreter','latex','FontSize',14);
set( findobj(gca,'type','line'), 'LineWidth', 2);
gsp_plotfig('Boundsparsity_path',paramplot2);



