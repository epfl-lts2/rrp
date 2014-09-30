%RR_INTRODUCTION

% This file produce the plot for the introduction
% No explanation is provided
% Author: Nathanael Perraudin



% clc;
clear;
close all;

global GLOBAL_save
global GLOBAL_baw

if GLOBAL_baw
    paramplot.pathfigure = 'figures/introduction/'; 
else
    paramplot.pathfigure = 'figures_color/introduction/';
end
paramplot.position = [100 100 300 200];
paramplot.titleweight ='bold';
paramplot.save=GLOBAL_save;
paramplot.baw=GLOBAL_baw; % plot everything in black and white

%% 


dynrange = 80;

% parameters
a=4;
M=1024;
Ltot=2*M;

Lwin = M;

% windows

% g = fir2long(fftshift(tukeywin(M/4,.4)),Lwin);
g = fir2long(gabwin('itersine',a,M/2),Lwin);

% ga = fir2long(gabwin('nuttall',a,M/8),Lwin);


figure;
plot_time(g);
title('(a)');
setplottime(g);
save_name='analysis_window';
plotfig(save_name,paramplot );
%%

%signal
f=sin(2*pi*((1:2048)/64).^2)';
f = f + sin(2*pi*(1:length(f))'/8);
f0 = f;
f(513:1536) = f(513:1536) + sin((pi*(1:1024)'/1024)).*sin(2*pi*(1:1024)'/4);

% spectrogramme
c=dgtreal(f,g,a,M);

figure;
plotdgtreal(c,a,M,1,dynrange);
title('(b)');
save_name='original_spectrogram';
plotfig(save_name,paramplot );

% spectrogramme
c0=dgtreal(f0,g,a,M);

figure;
plotdgtreal(c0,a,M,1,dynrange);
title('(c)');
save_name='denoised_spectrogram';
plotfig(save_name,paramplot );



%% Gabor synthesis
gd=gabdual(g,a,M,Lwin);
figure;
plot_time(gd);
title('(e)');
save_name='synthesis_windows_1';
plotfig(save_name,paramplot );
fprintf('  Reconstruction error of the canonical dual: %g \n',...
                                gabdualnorm(g,real(gd),a,M,Ltot));

%%

gd2 = real(gabfirdual(Lwin*2,g,a,M,'gamma',0,'mu',0,'var2t',1,...
        'maxit',100,'evolution',1));
% glike=ones(Lwin,1);
% gd2=real(gabglike(g,glike,a,M,'evolution',1,'support',Lwin));
%%
figure;
plot_time(long2fir(gd2,Lwin/4));
title('(g)');

save_name='synthesis_windows_2';
plotfig(save_name,paramplot );
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd2),a,M,Ltot));
                            

%%
cm=c;
cm(241:272,161:352)=0;
% 
% Weights = repmat(cos(2*pi*(1:501)/501),183, 1);
% cm(75:end,500:1000)= Weights.*cm(75:end,500:1000);

figure;
plotdgtreal(cm,a,M,1,dynrange);
title('(d)');
save_name='modified_spectrogram';
plotfig(save_name,paramplot );
%%

c1=dgtreal(idgtreal(cm,gd,a,M,Ltot),g,a,M);
c2=dgtreal(idgtreal(cm,gd2,a,M,Ltot),g,a,M);

% c1=dgtreal(idgtreal(cm,gd,a,M,Ltot),ga,a,M);
% c2=dgtreal(idgtreal(cm,gd2,a,M,Ltot),ga,a,M);

figure;
plotdgtreal(c1,a,M,1,dynrange);
title('(f)');
save_name='new_spectrogram_1';
plotfig(save_name,paramplot );
figure;
plotdgtreal(c2,a,M,1,dynrange);
title('(h)');
save_name='new_spectrogram_2';
plotfig(save_name,paramplot );


norm(c1(:)-c0(:))
norm(c2(:)-c0(:))

snr(c0,c1)
snr(c0,c2)

