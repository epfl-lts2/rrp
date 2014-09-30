function plot_ambiguity(g,dr)

z = 2.5;

L = length(g);
g = fir2long(g,ceil(z*L));
L = length(g);
G = abs(dgt(g,g,1,L));
G = G/max(G(:));
% norm(G(:),1);
G =20*log10(G);
%mG = max(G(:));

x=linspace(-z*L/2,z*L/2,z*L);
y=linspace(-z/2,z/2,z*L);

imagesc(x,y,fftshift(G),[-dr,0]);

colorbar;

%axis off;

end