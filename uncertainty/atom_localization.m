%ATOM_LOCALIZATION Illustration of the localization of some atoms
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
%   We present in Figure 1 the average and maximum hop distance
%   $h_{G}(i,\tilde{i})$. In this example, we control the concentration of
%   a kernel $\hat{g}$ with a dilation parameter $a$:
%   $\widehat{g_a}(x)=\hat{g}(ax)$. Increasing $a$ compresses the kernel in
%   the Fourier domain and increases the spread of the localized atoms in
%   the vertex domain. Note that even for high spectral compression, the
%   hop distance $h_{G}(i,\tilde{i})$ remains low. Additionally, we also
%   compute the mean relative error between $\|T_i g^2 \|_\infty$ and
%   $|T_ig^2(i)|$. This quantity asserts how well $\|T_i g \|_2^2$
%   estimates $\|T_i g^2\|_\infty$.
%
%   .. figure::
%
%      Localization experiment 
%
%      The heat kernel is defined as $\hat{g}(ax) = e^{-\frac{10\cdot ax}{\lambda_{\max}}}$
%      and the wavelet kernel $\hat{g}(ax) = \sqrt{40} \cdot ax \cdot e^{- \frac{40 \cdot ax}{\lambda_{\max} }}$. 
%      For a smooth kernel $\hat{g}$, the hop distance $h_{G}$ between $i$ 
%      and $\tilde{i}=\arg\max_j |T_ig(j)|$ is small.
%    
%   
%   References: perraudin2016global
%

clear
close all;

%% Plotparameter
global SAVE ;

paramplot.save = SAVE;
paramplot.position = [100 100 300 200];
paramplot.savefig = 0;


gsp_reset_seed(0)
%%


a = [0.1 0.2 0.5 1 2 5 10];
Nf = length(a);
N = 100;

G = gsp_sensor(N);
G = gsp_compute_fourier_basis(G);

g0 = gsp_design_heat(G,20);
g = cell(Nf,1);
for ii = 1:Nf
    g{ii} = @(x) g0(a(ii)*x);
end


atoms = gsp_vec2mat(gsp_filter_analysis(G,g,eye(N)),Nf);

ind = zeros(Nf,N);
err_ind = zeros(Nf,1);
err_m = zeros(Nf,1);
hope = zeros(Nf,1);
maxhope = zeros(Nf,1);
%%
for ii = 1:Nf
    [val,ind(ii,:)] = max(reshape(atoms(:,ii,:),N,N));
    for jj = 1:N
        tmp = gsp_hop_distanz(G,jj,ind(ii,jj));
        hope(ii) = hope(ii) + tmp/N;
        if maxhope(ii)<tmp;
            maxhope(ii) = tmp;
        end
    end
    err_ind(ii) = sum(abs(ind(ii,:)-(1:N))>0)/N;
    err_m(ii) = sum(abs((diag(reshape(atoms(:,ii,:),N,N))-val')./val'))/N*100;
end


%%

figure()
paramplot.show_sum = 0;
gsp_plot_filter(G,g,paramplot);
legend(num2str(a',2))
gsp_plotfig('scale_filter_heat',paramplot)

mat2tex([a',err_m,hope,maxhope],0,'%g')




