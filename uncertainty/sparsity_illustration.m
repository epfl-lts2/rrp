%SPARSITY_ILLUSTRATION Compute the sparsity level of some signals
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
%   Figure 1: uses some basic signals to illustrate this notion of
%   concentration, for different values of $p$.
%
%   .. figure::
%
%      Typical examples
%
%      
%      
%   
%   References: perraudin2016global
%

%%
close all;
clear


set(0,'DefaultTextInterpreter','Latex')

global SAVE ;

paramplot.save = SAVE;
paramplot.position = [100 100 300 220];
paramplot.savefig = 0;

%%
N=10;
f1=ones(N,1);
f1=f1/norm(f1);
f4=zeros(N,1);
f4(5)=1;
f3=zeros(N,1);
f3([2,3,8,9])=1/2;
G=gsp_path(10);
G=gsp_compute_fourier_basis(G);
EL=G.U*diag(exp(-G.e))*G.U';
f2=EL(:,5);
f2=f2/norm(f2);
F=[f1,f2,f3,f4];

figure(1);

for ii=1:4
    subplot(2,2,ii)
    stem(F(:,ii),'LineWidth',3);
    xlim([1,N]);
    ylim([0,1]);
%     set(gca,'FontSize',24);
    title(['$f_',num2str(ii),'$'])
end
gsp_plotfig('concentration_example',paramplot)

bounds=[1/sqrt(N),1]

s=@(f,p) (norm(f,p)/norm(f,2))*(p>2)+(norm(f,2)/norm(f,p))*(p<=2);



dat=zeros(5,5);
dat(:,1)=[1;4/3;2;4;Inf];
for ii=1:5
    for jj=2:5
        dat(ii,jj)=s(F(:,jj-1),dat(ii,1));
    end
end
% disp(dat)
fprintf('       p      s_p(f1)   s_p(f2)   s_p(f3)   s_p(f4)   \n');
mat2tex(dat,'table','%5.2f')