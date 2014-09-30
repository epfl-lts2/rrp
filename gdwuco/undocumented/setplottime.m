function [  ] = setplottime(g)
%SETPLOTTIME Set some parameters for frequency plot


max_g=max(g(:));
min_g=min(g(:));
range=max_g-min_g;
set(gca,'Ylim',[min_g-0.05*range max_g+0.05*range]);


end

