% This file produce the plot for the introduction
% No explanation is provided
% Author: Nathanael Perraudin

clc;
clear all;
close all;

paramplot.pathfigure = 'figures/presentation/'; 
paramplot.position = [100 100 300 200];
paramplot.titleweight ='bold';
paramplot.save=1;
paramplot.baw=1; % plot everything in black and white

%% Gabor transform

% parameters
a=8;
M=256;
Ltot=1024;

% windows
g=gabwin('nuttall',a,M);
g=fir2long(g,Ltot);



figure;
plot(fftshift(g));
title('Analysis window');
setplottime(Ltot,g);
save_name='analysis_window';
plotfig(save_name,paramplot );
%%

%signal
f=[ sin(2*pi*(1:256)/32),...
    sin(2*pi*(1:384)/128),sin(2*pi*((1:384)/100).^2)]';


figure;
plot(f)
title('Signal');
xlabel('Time')
setplottime(Ltot,f);
save_name='signal';
plotfig(save_name,paramplot );
%%

% spectrogramme
c=dgtreal(f,g,a,M);

figure;
plotdgtreal(c,a,M,1,50,'nocolorbar');
title('Spectrogram');
save_name='original_spectrogram';
plotfig(save_name,paramplot );


%% Gabor synthesis
gd=gabdual(g,a,M,Ltot);
figure;
plot(fftshift(gd));
title('Synthesis window');
setplottime(Ltot,gd);
save_name='synthesis_windows_1';
plotfig(save_name,paramplot );
%%
atemp=128;
gtemp=gabwin('nuttall',atemp,M);
gtemp=fir2long(gtemp,Ltot);
glike=gabdual(gtemp,atemp,M,Ltot);
glike=ones(512,1);
gd2=real(gaboptdual(g,a,M,'glike',glike,'omega',1));


figure;
plot(fftshift(gd2));
title('Synthesis window 2');
setplottime(Ltot,gd2);
save_name='synthesis_windows_2';
plotfig(save_name,paramplot );
%%
cm=c;
cm(62:66,62:66)=max(abs(c(:)));

figure;
plotdgtreal(cm,a,M,1,50,'nocolorbar');
title('Modified spectrogram');
save_name='modified_spectrogram';
plotfig(save_name,paramplot );
%%

c1=dgtreal(idgtreal(cm,gd,a,M,Ltot),g,a,M);
c2=dgtreal(idgtreal(cm,gd2,a,M,Ltot),g,a,M);

figure;
plotdgtreal(c1,a,M,1,50,'nocolorbar');
title('Spectrogram: modified signal 1');
save_name='new_spectrogram_1';
plotfig(save_name,paramplot );
figure;
plotdgtreal(c2,a,M,1,50,'nocolorbar');
title('Spectrogram: modified signal 2');
save_name='new_spectrogram_2';
plotfig(save_name,paramplot );

