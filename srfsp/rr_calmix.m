%RR_CALMIX Superresolution used in spectrometry for the data calmix
%   
%   Demonstration file on synthetic dataset
%   
%   Author:  Perraudin Nathanael
%
%   EPFL - LTS2  --  Mai 2014
%   
%   Dependencies
%   ------------
%
%   In order to use this matlab file you need the UNLocbox toolbox. You
%   can download it on https://lts2research.epfl.ch/unlocbox .
%   
%   The problem
%   -----------
%
%   In this file, we attempt to increase the resolution of a signal in the
%   Fourier domain. This signal is particular because we know that it is
%   sparse. The original signal in the spectral domain is displayed in the
%   figure 1.
%
%   .. figure::
%
%      Fig 1: Original signal in the spectral domain
%
%      We consider this signal as the ground true for our problem.
%
%   To obtain such a precision in frequency, we do have to capture the
%   signal for a long time. Since I do not have this amout of time to
%   finish my PHD thesis, I will reduce the time aquisition by a factor 20.
%   This is equivalent to a convolution in frequency with a sinc function.
%   The figures 2 and 3 provide a zoom in of the signal to demonstrate
%   this effect:
%
%   .. figure::
%
%      Fig 2: Original signal in the spectral domain - zoom in
%
%      We still have sharp pick in the Fourier domain
%
%   .. figure::
%
%      Fig 3: Measurements in the spectral domain - zoom in
%
%      Here we observe the effect of the convolution, the peaks are blured.
%
%   To recover the original signal, we will solve the following convex
%   optimization problem. More information can be found in
%   |rr_superresolution_for_spectrometry| :
%
%   .. argmin_x || x ||_1    such that || M F^{-1} x - y ||_2 < epsilon
%
%   .. math:: \operatorname{argmin}_x \|x\|_1   \text{ s. t. } \|M F^{-1} x - y \|_2 \leq \epsilon
%
%   where $y$ are the measurements, $F^{-1}$ the inverse Fourier transform
%   and $M$ the masking operator. $\epsilon$ is the radius of the
%   $l_2$-ball that is chosen with respect of the noise.
%
%   Solving the problem lead to the following deblured signal.
%
%   .. figure::
%
%      Fig 4: Deblured signal - zoom in
%
%      Solution of the convex optimization problem
%
%   We obtain something much more close to the original signal. At least,
%   the peak a much more sharp. For the amplitude some more work needs to
%   be done. We provide 3 more figures with a closer zoom in (Figures
%   5, 6, 7):
%
%   .. figure::
%
%      Fig 5: Original signal in the spectral domain - zoom in 2
%
%      Here we observe the spread of the peak of the original signal
%
%   .. figure::
%
%      Fig 6: Measurements in the spectral domain - zoom in 2
%
%      Here we observe the effect of the convolution, the peaks are not
%      peaks anymore
%
%   .. figure::
%
%      Fig 6: Deblured in the spectral domain - zoom in 2
%
%      We have peaks with comparable width as the original signal.
%


%% Initialisation

clear;
close all;

% Lauch the unlocbox
init_unlocbox;

% Parameters for the problem
Mult = 0.05;     % Measurements percents
filename = '2-calmix.hdf5'; % filename

verbose = 1;    % Verbosity level
maxit = 100;   % Maximum number of iteration
tol   = 10e-10;   % Tolerance to stop iterating

show_evolution = 0; % show evolution of the algorithm

epsilon = 0;% 1.1*sqrt(Ntot)*sigma;% Radius of the B2-ball


%% Load the data



info = h5info(filename);
%h5disp(filename);
signal = h5read(filename,['/',info.Datasets(2).Name]);
fs = h5read(filename,['/',info.Datasets(1).Name]);

Ntot = length(signal);

N = round(Ntot*Mult);

% Masking operation
Mask = zeros(Ntot,1);
Mask(1:N) = 1;
M = @(x) Mask.*x;

% mesurements
y = M(signal);
ffty = fftreal(y);

% Ground true
ffts = fftreal(signal);





%% Setting the problem in the unlocbox

% data fidelity term

A = @(x) fftreal(M(x));
At = @(x) M(ifftreal(x,Ntot));
paramfid.A = At;
paramfid.At = A;
paramfid.nu = 1;
paramfid.tight = 1;
paramfid.verbose = verbose - 1;
paramfid.y = y;
paramfid.epsilon = epsilon;

ffid.eval = @(x) eps;
ffid.prox = @(x,T) proj_b2(x,T,paramfid);

% Prior term
paraml1.verbose = verbose - 1;
fl1.eval = @(x) norm(x,1);
fl1.prox = @(x,T) prox_l1(x,T,paraml1);


%% Solving the problem
paramsolver.gamma = max(abs(ffty));
paramsolver.maxit = maxit;
paramsolver.tol = tol;
paramsolver.verbose = verbose;

if show_evolution
    fig=figure(100);
    paramsolver.do_sol=@(x) plot_signal_fftreal(x,10,fig,fs,ffts);
end

%x0 = guess_sol(ffty);
x0 = zeros(size(ffty));
sol1 = douglas_rachford(x0,fl1,ffid,paramsolver);




guess = guess_sol(sol1);
sol = douglas_rachford(guess,fl1,ffid,paramsolver);


if show_evolution
    close(fig)
end



%% Display the results

figure;
plotfftreal(ffts,fs);
title('Ground true');


% figure;
% plotfftreal(ffty,fs);
% title('Mesurements');
% 
% figure;
% plotfftreal(sol,fs);
% title('Recovered signal');



figure;
plotfftreal(ffts,fs);
xlim([1.6e5,3e5]);
title('Ground true zoom in');


figure;
plotfftreal(ffty,fs);
xlim([1.6e5,3e5]);
title('Mesurements zoom in');

figure;
plotfftreal(sol,fs);
xlim([1.6e5,3e5]);
title('Recovered signal zoom in');


figure;
plotfftreal(ffts,fs);
axis([1.643e5,1.649e5,0,3e7]);
title('Ground true zoom in 2');


figure;
plotfftreal(ffty,fs);
axis([1.643e5,1.649e5,0,3e7]);
title('Mesurements zoom in 2');

figure;
plotfftreal(sol,fs);
axis([1.643e5,1.649e5,0,3e7]);
title('Recovered signal zoom in 2');




