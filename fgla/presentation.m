% This file produce the file for the WASPAA presentation in NY the 23rd
% October 2013
%
% Author: Nathanael Perraudin
% Date  : 11st October 2013

%% Initialisation
clear;
close all;

paramplot.pathfigure = 'figures/presentation/'; 
paramplot.position = [100 100 300 200];
paramplot.titleweight ='bold';
paramplot.save=0;
paramplot.baw=0; % plot everything in black and white

%% Gabor transform

%% Gabor transform

% parameters
a=1;
Ltot=256;
M=Ltot;


% windows
g=gabwin('nuttall',a,64);
g=fir2long(g,Ltot);
gd=gabdual(g,a,M);


figure;
plot(fftshift(g));
title('Analysis window');
setplottime(g);
save_name='analysis_window';
plotfig(save_name,paramplot );


%signal
t=1:Ltot;
f=sin(2*pi*t.^3/Ltot^2/3);
maskt=zeros(size(f));
gt=gabwin('itersine',a,100);
maskt(41:140)=5*fftshift(gt);
f=f+maskt.*sin(2*pi*t*0.4);


figure;
plot(f)
title('Signal');
xlabel('Time')
setplottime(f);
save_name='signal';
plotfig(save_name,paramplot );


% spectrogramme
c=dgtreal(f,g,a,M);

figure;
plotdgtreal(c,a,M,1,50,'nocolorbar');
title('Spectrogram');
save_name='original_spectrogram';
plotfig(save_name,paramplot );


% new spectrogram

Mask=ones(size(c));
Mask(85:120,26:135)=0;
Mask(85:120,190:210)=0;
Mask(10:14,:)=0;
c2=c.*Mask;
figure;
plotdgtreal(c2,a,M,1,50,'nocolorbar');
title('Modified spectrogram');
save_name='spectrogram_hole';
plotfig(save_name,paramplot );
figure;
plotdgtreal(Mask,a,M,1,50,'nocolorbar');
title('Mask');
save_name='mask';
plotfig(save_name,paramplot );

% Synthesis 
f2=idgtreal(c2,gd,a,M);
c3=dgtreal(f2,g,a,M);
figure;
plotdgtreal(c3,a,M,1,50,'nocolorbar');
title('New spectrogram');
save_name='new_spectrogram_hole';
plotfig(save_name,paramplot );
figure;
plot(f2)
title('New signal');
xlabel('Time')
setplottime(f2);
save_name='signal_filtered';
plotfig(save_name,paramplot );
% remove the phase
c4=abs(c);
% Synthesis 
f4=idgtreal(c4,gd,a,M);
c5=dgtreal(f4,g,a,M);
figure;
plotdgtreal(c5,a,M,1,50,'nocolorbar');
title('New Spectrogram');
save_name='new_spectrogram_nophase';
plotfig(save_name,paramplot );
figure;
plot(f4)
title('Reconstructed signal');
xlabel('Time')
setplottime(f4);
save_name='signal_reconstructed';
plotfig(save_name,paramplot );

%% Plot some other spectrograms
sig=bat;
figure;
sgram(sig);
plotfig('bat',paramplot);

sig=gspi;
figure;
sgram(sig);
plotfig('gspi',paramplot);


break;
%%
% Use GLA
[~,GB]= gabframebounds(g,a,M);
G = @(x) dgtreal(x,g,a,M);
Gt= @(x) idgtreal(x,g,a,M)/GB;


param.maxit=100;
param.verbose=1;
param.tol=10e-4;
alpha=0.99;
[ f5,c6,info_reconstruct1] = spectrogram_reconstruction( c4,G,Gt,...
                                     'GLA',alpha, param);
figure;
plotdgtreal(c6,a,M,1,50,'nocolorbar');
title('Griffin Lim Algorithm');
save_name='new_spectrogram_nophase_gla';
plotfig(save_name,paramplot );                                

                                 
                                 
[ f6,c7,info_reconstruct2] = spectrogram_reconstruction( c4,G,Gt,...
 'FGLA',alpha, param);
figure;
plotdgtreal(c7,a,M,1,50,'nocolorbar');
title('Fast Griffin Lim algorithm');
save_name='new_spectrogram_nophase_fgla';
plotfig(save_name,paramplot );

figure
R(1,:)=info_reconstruct1.ssnr;
R(2,:)=info_reconstruct2.ssnr;
iterations=1:param.maxit;
semilogx(iterations,R);
legend('Griffin Lim','Fast Griffin Lim','Location','SouthEast');
xlabel('Iterations ','FontSize',12);
ylabel('ssnr (dB)','FontSize',12);
save_name='convergence_comparaison';
plotfig(save_name,paramplot );
drawnow;
