%% Demo introduction

clear;
close all;

gsp_reset_seed(3);

%% create a band limited stationary signal

fs = 1;         % sampling frequency
tmax = 50;      % maximum time
fc = 0.05;      % center frequency for the frequency band

N = tmax*fs;    % length of the signal

% frequency axis
fa = ifftshift(linspace(-fs/2,fs/2-fs/N,N));
% time axis
t = 1/fs:1/fs:tmax;



% band center frequency 
sigma = 0.02;
bf = @(x) exp(-(x-fc).^2/sigma.^2) + exp(-(x+fc).^2/sigma.^2);
filt = bf(fa);
f = @(x) real(ifft(fft(x).*filt));
% 
% figure(4)
% stem(fa,filt)
% xlabel('Normalized frequency')
% ylabel('Amplitude')

% Generate the signal
w = randn(1,N);
s = f(w);


%% Find another possible end for the signal

tlim = 43;
Nlim = tlim *fs;


sigma = 0.02;
bf2 = @(x) exp(-(x-3*fc).^2/sigma.^2) + exp(-(x+3*fc).^2/sigma.^2);
filt2 = bf2(fa);
f2 = @(x) real(ifft(fft(x).*filt2));

mask = zeros(1,N);
mask( 1:Nlim ) = 1;
M = @(x) mask.*x;
obs = M(s);

sr = randn(1,N);
for ii = 1:300
    sr = f2(sr);
    sr(logical(mask)) = obs(logical(mask));
end


sunlike = sr;
%sunlike(logical(1-mask)) = -sr(logical(1-mask))+2*sr(Nlim);

t1 = t(1:Nlim);
s1 = s(1:Nlim);

figure(1)
plot(t,sunlike,'go--',t,s,'ro-',t1,s1,'bo-');
legend('Prediction 1', 'Prediction 2','Observed value','Location','SouthWest')
xlabel('Time (s)')
ylabel('Signal')
paramplot.position = [100 100 500 280];
gsp_plotfig('intro_stationary_signal',paramplot)


%%
gsp_reset_seed(6)
N = 100;
G = gsp_spiral(N,2);
G = gsp_compute_fourier_basis(G);

fc = G.lmax/20;
sigma = G.lmax/10;
g = @(x) exp(-(x-fc).^2/sigma.^2).*exp(-x/sigma*5);

w = randn(N,4);
s = gsp_filter(G,g,w);

figure(2)

paramplot.show_edges = 0;

subplot(221)
title('Signal 1')
gsp_plot_signal(G,s(:,1),paramplot)
colorbar off
subplot(222)
title('Signal 2')
gsp_plot_signal(G,s(:,2),paramplot)
colorbar off
subplot(223)
title('Signal 3')
gsp_plot_signal(G,s(:,3),paramplot)
colorbar off
subplot(224)
title('Signal 4')
gsp_plot_signal(G,s(:,4),paramplot)
colorbar off


paramplot.position = [100 100 300 280];
%title('Signal on an 2d irregular domain')
gsp_plotfig('intro_irr_signal',paramplot)

figure(3)
gsp_plot_signal(G,s(:,1))
title('Signal 1 on the Graph')
gsp_plotfig('intro_graph',paramplot)


