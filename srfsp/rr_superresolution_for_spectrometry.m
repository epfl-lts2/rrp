%RR_SUPERRESOLUTION_FOR_SPECTROMETRY Superresolution used in spectrometry
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
%   This last signal pocess a very high resolution frequency resolution. To
%   obtain such a frequency resolution, a very long time measurement is
%   done. Let's suppose now, that we would like to reduce this measurment
%   time to something shorter. This is equivalent to convolute the signal
%   in the Fourier domain by a sinc. Figure 2 show the effect of this
%   convolution. This operation is considered as masking of all unmeasured
%   coefficients.
%
%   .. figure::
%
%      Fig 2: Blured signal
%
%      This signal is the result of the convolution with a sinc 
%
%   Finally, we add some gaussian noise to test if our method is robust.
%   Figure 3 display the noiy measurments that we use for our problem.
%
%   .. figure::
%
%      Fig 3: Noisy blured signal
%
%      Measurement used for our problem
%
%   To get close to the ground true, we will solve the following convex
%   optimization problem:
%
%   .. argmin_x || F x ||_1    such that || M x - y ||_2 < epsilon
%
%   .. math:: \operatorname{argmin}_x \|F x\|_1   \text{ s. t. } \|M  x - y \|_2 \leq \epsilon
%
%   where $y$ are the measurements, $F^{-1}$ the inverse Fourier transform
%   and $M$ the masking operator. $\epsilon$ is the radius of the
%   $l_2$-ball that is chosen with respect of the noise.
%
%   Solving the problem lead to the following denoised signal.
%
%   .. figure::
%
%      Fig 4: Denoised signal
%
%      Solution of the convex optimization problem
%
%   The problem is solved in two times. We compute a first estimate by
%   solving the convex problem. Then with this rougth estimate, we guess
%   the position and the amplitude of the dirac to obtain a new estimation.
%   This estimation is then reinsterted as starting point into the convex
%   problem to obtain a final solution.
%


%% Initialisation

clear;
close all;

% Lauch the unlocbox
init_unlocbox;
ltfatstart;

% Parameters for the problem
Ns = 5;         % Number of sines
Amp = [1,2];  % Minimum/Maximum amplitude for the sine
fs = 1000;      % Sampling frequency
T = 1;          % Sampling time
Ttot = 100;     % Total time
sigma = 0.5;      % Noise level

Ntot = fs*Ttot;

epsilon = 1.1*sqrt(fs*T)*sigma;% Radius of the B2-ball

verbose = 1;    % Verbosity level
maxit = 1000;   % Maximum number of iteration
tol   = 10e-10;   % Tolerance to stop iterating

show_evolution = 1; % show evolution of the algorithm


%% Creation of the problem


s = zeros(Ntot,1);

% create the sines
for ii = 1:Ns
    fii = round(rand(1)*Ntot) / Ntot;
    ampii = Amp(1) + rand(1) * (Amp(2)-Amp(1));
    s = s + +ampii * sin( 2 * pi * fii * (1:Ntot)' );
    
end

% add the noise
sn = s + sigma*randn(Ntot,1);

% Masking operation
Mask = zeros(Ntot,1);
Mask(1:fs*T) = 1;
M = @(x) Mask.*x;

% mesurements
y = M(sn);
ffty = fftreal(y);

% Ground true
ffts = fftreal(s);



%% Setting the problem in the unlocbox

% data fidelity term


paramfid.A =  M;
paramfid.At = M;
paramfid.nu = 1;
paramfid.tight = 1;
paramfid.verbose = verbose - 1;
paramfid.y = y;
paramfid.epsilon = epsilon;

ffid.eval = @(x) eps;
ffid.prox = @(x,T) proj_b2(x,T,paramfid);

% Prior term
A = @(x) fftreal(x)/sqrt(Ntot);
At = @(x) ifftreal(x,Ntot)*sqrt(Ntot);
paraml1.A = A;
paraml1.At = At;
paraml1.tight = 1;
paraml1.verbose = verbose - 1;
fl1.eval = @(x) norm(A(x),1);
fl1.prox = @(x,T) prox_l1(x,T,paraml1);


%% Solving the problem
paramsolver.gamma = 1 ;
paramsolver.maxit = maxit;
paramsolver.tol = tol;
paramsolver.verbose = verbose;

if show_evolution
    fig=figure(100);
    paramsolver.do_sol=@(x) plot_signal_fftreal_2(x,10,fig,fs,ffts);
end

sol = douglas_rachford(y,fl1,ffid,paramsolver);
%guess = guess_sol(sol1);

%%
%sol = douglas_rachford(guess,fl1,ffid,paramsolver);





%% Improve the solution with l12

% % Prior term
% tg = 20;
% tau = 0.00001;
% paraml12.verbose = verbose - 1;
% fl12.eval = @(x) tau*norm_l12(reshape(x(1:end-1),tg,[])');
% fl12.prox = @(x,T) [reshape(transpose(prox_l12(transpose(reshape(x(1:end-1),tg,[])),tau*T,paraml12)),[],1);0];
%
% sol = ppxa(sol,{fl1,ffid,fl12},paramsolver);



%% Display the results

figure;
plotfftreal(ffts,fs);
title('Ground true');

fftsm = fftreal(M(s));
figure;
plotfftreal(fftsm,fs);
title('Mesurements without noise');

figure;
plotfftreal(ffty,fs);
title('Mesurements');

figure;
plotfftreal(A(sol),fs);
title('Recovered signal');




