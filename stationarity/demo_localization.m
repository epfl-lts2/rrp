%% Illustration of the localization operator

clear
close all

gsp_reset_seed(0);
N = 100;
G = gsp_sensor(N);
G = gsp_compute_fourier_basis(G);

nodes = [1,34,83];

g = @(x) (10*x)/G.lmax.*exp(-(10*x).^2/G.lmax^2);



figure(1)
subplot(221)
gsp_plot_filter(G,g);
title('Filter, Mexican hat')

c = [-0.2 0.5];
subplot(222)
paramplot.vertex_highlight = nodes(1);
gsp_plot_signal(G,gsp_localize(G,g,nodes(1)),paramplot);
caxis(c);
subplot(223)
paramplot.vertex_highlight = nodes(2);
gsp_plot_signal(G,gsp_localize(G,g,nodes(2)),paramplot);
caxis(c);

subplot(224)
paramplot.vertex_highlight = nodes(3);
gsp_plot_signal(G,gsp_localize(G,g,nodes(3)),paramplot);
caxis(c);

gsp_plotfig('demo_localization')