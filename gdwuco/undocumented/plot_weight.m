%PLOT_WEIGHT


%% initialization
%clear all;
close all;


global GLOBAL_save
global GLOBAL_baw
 

% plotting
dr=100;
if GLOBAL_baw
    paramplot.pathfigure = 'figures/'; 
else
    paramplot.pathfigure = 'figures_color/'; 
end
    
paramplot.position = [100 100 300 200];
paramplot.titleweight ='bold';
paramplot.save = GLOBAL_save; % save figures
paramplot.baw = GLOBAL_baw; % plot everything in black and white


%%
L = 240;


if mod(L,2)
    w=[0:1:(L-1)/2,(L-1)/2:-1:1]';
else         
    w=[0:1:L/2-1,L/2:-1:1]';
end
w=w/sqrt(L);
W = repmat(w,1,L).^2 + repmat(w',L,1).^2;   
W=log(1+W);

x=linspace(-L/2,L/2,L);
y=linspace(-1/2,1/2,L);


imagesc(x,y,fftshift(W))
colorbar
title('Weights')
xlabel('Time')
ylabel('Frequency')
plotfig('weights',paramplot)


