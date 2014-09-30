%RR_GABOPTDUAL Optimization of dual gabor windows
%                           
%   Reproducible research addendum for optimization of dual gabor windows
%   ---------------------------------------------------------------------
%   
%   DESINGING GABOR WINDOWS USING CONVEX OPTIMIZATION
%
%   Paper: Nathanael Perraudin, Nicki Holighaus, Peter L. Sondergaard and
%   Peter Balazs 
%   
%   Demonstration matlab file:  Perraudin Nathanael
%
%   ARI -- January 2013
%   
%   Dependencies
%   ------------
%
%   In order to use this matlab file you need the UNLocbox toolbox. You
%   can download it on http://unlocbox.sourceforge.net and the LTFAT
%   toolbox. You can download it on http://ltfat.sourceforge.net
%
%   You can also want the convergence of the algorithm at
%   http://unlocbox.sourceforge.net/rr/gdwuco/gif.php.
%
%   Goal
%   ----
%
%   We show the deference of different dual sythesis windows designed
%   either by classic techniques or by convex optimization. 
%
%   All the examples provided in this experiment were computed with an
%   Itersine analysis window with $L_g=60$ $L=240$, $a=15$ and $M=120$,
%   without support constraints. This setup, in particular its high
%   redundancy, allows us to shape the dual windows rather freely for
%   different objective functions, therefore producing characteristic
%   examples. The window is shown below
%
%   .. figure::
%
%      Analysis window in time domain
%
%      The chosen windows is an itersine of length 60
%
%   .. figure::
%
%      Analysis window in frequency domain
%
%      
%
%   The canonical dual, equals the window up to scaling.  
%
%   .. figure::
%
%      Canonical dual window in time domain
%
%       
%
%   .. figure::
%
%      Canonical dual window in frequency domain
%
%      
%
%   To obtain different dual window than the canonical dual, we propose to
%   solve a convex optimization problem of the form:
%
%   .. gd  = argmin_x    f(x)
%
%   ..     such that    x is a dual windows of g and is compaclty supported
%
%   .. math:: \begin{split}  \text{gd}  =  \text{arg} \min_x   & f(x) \\     \text{such that }& x \text{ is a dual windows of }g  \end{split}
%
%   with 
%         - $x$         :  optimization variable
%         - $f(x)$      :  functional to promote properties.
%
%   We compare several functions to be minimized. A list is presented here
%   below. To obtained the best window, we may use several of them:
%   
%   * $|| \nabla x ||_2^2$  : force the windows to be smooth
%   
%   * $|| \nabla F x ||_2^2$: force the windows to be localized
%
%   * $|| x ||_2^2$         : minimize the two norm of the synthesis windows
%
%   * $|| x ||_1$           : minimize the one norm of the synthesis windows
%
%   * $|| F x ||_1$         : minimize the one norm in frequencie of the windows
%
%   * $|| x - z ||_2^2$     : Make x look like z in the 2 norm sense
%
%   * $var(x)$              : Localize x in time
%
%   * $var(x^2)$            : Localize x in time
%
%   * $var(F x)$            : Localize x in frequency
%
%   * $var((F x)^2)$        : Localize x in frequency
%
%   Solving the problem for each of these prior at once allows to produce
%   characteristic example to illustrate our technique
%
%   Lp-norm priors
%   --------------
%
%   In order to tune the solution of a convex optimization problem towards
%   the properties we desire, we have to select priors $f$ that
%   promote these properties. In this contribution, we mostly consider
%   priors that are fairly standard in optimization, or simple extensions
%   of such priors. Although their various effects are quite well known,
%   e.g. $l_1$ optimization favors solutions with a few large values
%   while an $l_2$ prior favors a more even spread of the energy,
%   considerable limitations are imposed by the duality constraint.
%   Although the set of Gabor dual windows is characterized by the WR
%   equations, their implications in terms of window shape, localization,
%   decay etc. remain largely unexplored. Therefore a short discussion of
%   relevant priors, expected effects and their actual effect in our
%   context seems worthwhile.  
%
%   An $l_2$ prior will, in our context, always lead to the canonical dual
%   Gabor window. In general, this prior will affect the values in a more
%   proportional way over the whole signal range. It is traditionally used
%   as a data fidelity term, i.e. the solution is expected to be close, in
%   the $l_2$-norm sense, to a given estimate. The associated objective
%   function is not only convex, but also smooth, admitting gradient
%   descent approaches for minimization.
%
%   .. figure::
%
%      Synthesis dual window in time domain minimizing the l2 norm
%
%      If we minimize the two norm of the sythesis window, we recover the
%      canonical dual. Here the original window.
%
%   .. figure::
%
%      Synthesis dual window in frequency domain minimizing the l2 norm
%
%      If we minimize the two norm of the sythesis window, we recover the
%      canonical dual.  
%
%   $l_1$-norm minimization is usually considered, whenever (approximate)
%   sparsity, i.e. a small number of (significant) non-zero values is
%   desired. As the convex relaxation of the $l_0$ minimization problem, it
%   is equivalent or at least close to sparsity optimization under certain
%   conditions. In general, these conditions are not satisfied by the WR
%   eq. Nevertheless, the restrictions imposed by the WR equations usually
%   allow a solution with few large values. It should be noted however,
%   that small $l_1$-norm does not imply clustering of the large values,
%   i.e. a solution supported on a short interval. As expected, some
%   concentration is induced by the duality constraint, see the figure
%   bellow. In the presented experiment, the $l_1$ solution possesses only
%   $15$ values above $-80$~dB (relative to the maximum amplitude) on an
%   interval around $0$, only half the number of WR equations $L a/M=30$.
%   However, other configurations have provided solutions with few
%   significant values spread over a larger interval. 
%
%   The proximity operator of the $\ell^1$ prior is computed by soft-thresholding:
%   
%   ..     soft_mu(y)  = sign(y) ( |y| - mu )_+
%   
%   .. math::  \text{soft}_\mu(y)=\text{sgn}(y) \left(|y|-\mu \right)_+ 
%
%   where $(.)_+=\max(.,0 )$. For compactly supported windows, a strictly
%   bandlimited dual window is usually not feasible. Therefore, when
%   applied in the Fourier domain, the $l_1$ prior cannot achieve a truly
%   sparse solution, but promotes a small number of significant values. In
%   many cases, the result is similar to actual concentration measures.
%
%   .. figure::
%
%      Synthesis dual window in time domain minimizing the l1 norm
%
%      This window minimize the one norm in time. It tend to have a few
%      big coefficient and all other equal to zero.
%     
%   .. figure::
%
%      Synthesis dual window in frequency domain minimizing the l1 norm
%
%      The window spread a lot in frequency. Indeed no optimization is done
%      in this domain.
%
%   .. figure::
%
%      Synthesis dual window in time domain minimizing the l1 norm
%
%      We provide a log scale figure.
%
%   .. figure::
%
%      Synthesis dual window in time domain minimizing the l1 norm in frequency
%
%      This looks very similar to a smoothness constraint. It is not
%      necessary the case. It is due to the WR equation system.
%
%   .. figure::
%
%      Synthesis dual window in frequency domain minimizing the l1 norm in frequency
%
%      The one norm in frequency is minimized. The aglorithm try to
%      sparsify the Fourier transform the window     
%
%   Concentration inducing function
%   -------------------------------
%
%   Our main objective in the following in this contribution will be the search for a
%   Gabor dual window with optimized/modified TF concentration. Therefore
%   we recall a number of different concentration measures. Inspired by the
%   famous Heisenberg inequality, the most natural way to impose
%   localization is to optimize the variance of the signal $x\in \mathbb{R}^L$ or
%   more precisely, its modulus:     
%
%   .. var(x) = 1/ sqrt(L)  sum_{i=-L/2}^{L/2-1} (i-\bar{|x|})^2 | x_i |
%
%   .. math: \text{var}(|x|)= 1/\sqrt{L} \sum_{i=-L/2}^{L/2-1}(i-\overline{|x|})^2|x_i|
%
%   with $\bar{|x|}=\sum_{i=-L/2}^{L/2-1}i | x_i |$ being the center of
%   gravity. As we consider symmetric windows $\bar{|x|}=0$, we can
%   simplify this expression to: $ var(|x|)= 1/\sqrt{L}
%   \sum_{i=-L/2}^{L/2-1}(i)^2 | x_i |$. In that case, the variance turns out
%   to be a weighted $l_1$-norm with quadratic weight $w^2$, 
%   $w := 1 / \sqrt{L} [-L/2, ..., L/2-1]$. Compared to $l_1$ minimization, 
%   this prior additionally penalizes values far from the origin, inducing
%   concentration. The proximity operator of $var(|x|)$ is a variation of
%   the $l_1$ proximity operator and computed by weighted soft
%   thresholding. And example is shown in the figure below.
%
%   .. figure::
%
%      Synthesis dual window in time domain minimizing the variance it time
%
%      We obtain almost a rectangle
%
%   .. figure::
%
%      Synthesis dual window in frequency domain minimizing the variance it time
%
%      
%
%   We also consider the variance of the energy of the signal:
%   $var(|x|^2)$, for symmetric windows equal to a weighted $l_2$ norm with
%   linear weight $w$: $var(x^2)=||w  x||_2^2$. Explicit computation of the
%   proximity operator leads to
%
%   .. prox_{gamma var(x^2) }(y) = 1 / ( 1 + 2 gamma w^2) * y
%   
%   .. math: \text{prox}_{\gamma \text{var}(x^2) }(y) = \frac{1}{1+2\gamma w^2} y,
%
%   i.e. multiplication with a function that decays quadratically away from
%   zero. The optimazation result is shown in the figure below.
%
%   .. figure::
%
%      Synthesis dual window in time domain minimizing the energy variance it time
%
%      There is also a good concentration in time, but with a kind of
%      smoothness.
%
%   .. figure::
%
%      Synthesis dual window in frequency domain minimizing the energy variance it time
%
%      In frequency, we have a better decay compared to the optimization of
%      the variance only.
%
%   A closely related concentration measure is smoothness in frequency, as
%   measured by the gradient of the Fourier transform $||\nabla F
%   x||_2^2$. Indeed, the resulting proximity operator has almost the same
%   form: 
%   
%   .. prox_{gamma ||\nabla F x||_2^2 }(y) = 1 / ( 1 + 2 gamma psi) * y
%
%   .. math \text{prox}_{\gamma \|\nabla \mathcal{F} x\|_2^2 }(y) = \frac{1}{1+2\gamma \psi} y
%
%   with $\psi[l] = 2-2 \cos(2 \pi l / L)$. Since $\psi[l] \approx C l^2$
%   for small $l$ and values away from $0$ are strongly attenuated, the
%   priors $var(|x|^2)$ and $||\nabla F x||_2^2$ often lead to similar
%   results. Both functions induce concentration by attenuation of values
%   far from the origin. Examples are shown in the figure below
%
%   .. figure::
%
%      Synthesis dual window in time domain minimizing the gradient of the Fourier transform
%
%      We observe that smoothing the frequencies result in concentrating the
%      window in time.
%
%   .. figure::
%
%      Synthesis dual window in frequency domain minimizing the gradient of the Fourier transform
%
%      
%
%   Concentration in frequency is easily achieved through $var(|| F x||)$,
%   $var(|| F x ||^2)$ or $||\nabla x||_2^2$. The respective proximity
%   operators are obtained simply by conjugating the proximitiy operators
%   discussed above with the (inverse) Fourier transform. The figure are
%   presented bellow.
%
%   .. figure::
%
%      Synthesis dual window in time domain minimizing the variance in frequency
%
%      
%
%   .. figure::
%
%      Synthesis dual window in frequency domain minimizing the variance in frequency
%
%      
%
%   .. figure::
%
%      Synthesis dual window in time domain minimizing the energy variance in frequency
%
%      This is equivalent to smooth the window in time.
%
%   .. figure::
%
%      Synthesis dual window in frequency domain minimizing the energy variance in frequency
%
%
%
%   .. figure::
%
%      Synthesis dual window in time domain minimizing the gradient in time
%
%      This windows is really smooth. In fact, it is the smoothest dual
%      window for the chosen parameter.
%
%   .. figure::
%
%      Synthesis dual window in frequency domain minimizing the gradient in time
%
%      The content in frequency of this windows is more concentrated in
%      low frequency. Indeed, optimizing smoothness concentrate the window
%      in low frequency by applying low pass filters.
%
%
%   Simulataneous time/frequency concentration function
%   ---------------------------------------------------
%
%   For simultaneous concentration in time and frequency, we can consider
%   jointly the time- and frequency-domain variants of the priors discussed
%   above. Alternatively, we use a single cost functions providing
%   concentration in both domains at once. In TF literature, modulation
%   space norms, i.e. $l_p$-norms on the short-time Fourier coefficients
%   are frequently used to measure joint TF localization.. In particular
%   $||x||_{S_0}=|| G_{g,1,L} x||_1$, where $g$ is a Gaussian function, is
%   considered as quality measure for window functions. Minimization of the
%   $S_0$-norm can be expected to yield TF concentrated windows. However,
%   similar to the $\ell^1$-norm, small $S_0$-norm does not guarantee
%   concentration around the origin (or any single TF location). As
%   expected, the duality constraint seems to enforce some localization,
%   though. 
%
%   Compared to the previously discussed priors, $S_0$ optimization is
%   considerably more expensive. Since we are not aware of a explicit
%   solution to the $S_0$ proximity operator, we propose its computation
%   via an iteration based on ADMM. The number
%   of required ADMM steps per PPXA iteration is low and scales well with
%   $L$ (usually 3-4 steps provided sufficient precision), but each substep
%   requires the computation of one full STFT and one inverse STFT, with a
%   complexity of $O(L^2 \log(L))$ each. An example for $S_0$
%   optimization is shown in the figure below.        
%
%   .. figure::
%
%      Synthesis dual window in time domain minimizing the S0 norm
%
%      This window is also well concentrated in time.
%
%   .. figure::
%
%      Synthesis dual window in frequency domain minimizing the S0 norm
%
%      A well concentrated windows in frenquency
%
%   In some cases, concentration can be further increased by $S_0$-norm
%   instead. The proximity operator is realized similar to the unweighted
%   case. The next figures show an example, considering the following
%   circular weights:
%   
%   ..          W[f,t]= ln( 1 + w^2[t] + w^2[f] )
%
%   .. math: W[f,t]=\ln \left( 1 + w^2[t]+w^2[f] \right),
%   
%   using the weight $w$ as defined above.
%
%   While other weights are clearly feasible, the weight above has been
%   tuned to yield good results in our experiments.
%
%   .. figure::
%
%      Synthesis dual window in time domain minimizing a weighted S0 norm
%
%      This window is also well concentrated in time.
%
%   .. figure::
%
%      Synthesis dual window in frequency domain minimizing a weighted S0 norm
%
%      A well concentrated windows in frenquency
%
%   Other cost functions
%   --------------------
%
%   The list of possible cost functions is vast and the full exploration of
%   the possibilities of convex optimization in window design is far beyond
%   the scope of any single contribution. As a rather academic example, we
%   propose a free design approach that selects the dual Gabor window
%   closest to the linear span of a model window $g_{sh}$, i.e. we find     
%
%   ..       argmin || x - P_<g_sh>(x)||_2    such that x is dual
%
%   .. math:  \mathop{\operatorname{arg~min}}\limits _{x \in \mathcal{C}_{\text{dual}}}\|x - P_{\langle g_{s}\rangle}x\|_2^2,
%
%   where $< g_s>$ is the linear span of $g_s$. The solution is computed by
%   a POCS (projection onto convex set) algorithm. Due to the examples
%   academic nature, we were not concerned with convergence time. Examples
%   using a sine wave and a dirac pulse as model window are presented in
%   the figure bellow
%
%   .. figure::
%
%      Objectiv shape function: a sine
%
%      We choose a sine function as an academic example.
%
%   .. figure::
%
%      Synthesis dual window in time domain
%
%      This window looks like a sine and it is still dual. 
%
%
%   .. figure::
%
%      Shape objectiv function: a dirac
%
%      We choose a dirac function. Why not?
%
%   .. figure::
%
%      Synthesis dual window in time domain
%
%      This window looks a bit like a dirac and it is still dual. 
%
%   .. figure::
%
%      Synthesis dual window in time domain (zoom in)
%
%      Note that the solution window (d) is actually composed of a smooth
%      bump function in addition to the 3 clearly visible impulses. 
%



%   Author : Nathanael Perraudin
%   Date   : Feb 24 2013




%% initialization
%clear all;
close all;


global GLOBAL_save
global GLOBAL_baw
global GLOBAL_evolution
 

% plotting
dr=100;
if GLOBAL_baw
    paramplot.pathfigure = 'figures/gaboptdual/'; 
else
    paramplot.pathfigure = 'figures_color/gaboptdual/'; 
end
    
paramplot.position = [100 100 300 200];
paramplot.titleweight ='bold';
paramplot.save = GLOBAL_save; % save figures
paramplot.baw = GLOBAL_baw; % plot everything in black and white

evo = GLOBAL_evolution;

%%

% general parameter
maxit=200; % maximum number of iteration
tol=10e-6; % tolerance to stop iterating

% parameters of frame and length
a=15;
M=120;
L=240;

% windows
g=firwin('itersine',60);
%g = fftshift(tukeywin(60,.6));

g=fir2long(g,L);



% original window

figure()
plot_time(g);
title('(a)');
save_name='analysis_time';
plotfig(save_name,paramplot );

figure();
magresp(g,'fir','dynrange',dr,'1');
title('(b)');
setplotfreq();
save_name='analysis_freq';
plotfig(save_name,paramplot );


% canonical dual
gcan=gabdual(g,a,M);
figure();
plot_time(gcan);
title('(c`)');
save_name='can_dual_time';
plotfig(save_name,paramplot );


figure()
magresp(gcan,'fir','dynrange',dr,'1');
title('(c`)');
setplotfreq();
save_name='can_dual_freq';
plotfig(save_name,paramplot );



%% Convex optimization L2 norm constraint

%parameter for convex optimization

omega=1; % weight for the L2 norm
mu=0;      % smoothing parameter
gamma=0;

drawnow;
% solving the problem
fprintf(' Solve the convex optimization problem: L2 norm\n');
gd=gaboptdual(g,a,M,'omega',omega,'mu',mu,...
    'gamma',gamma,'maxit',maxit,'tol',tol,'evolution',evo);
gd=real(gd);

[crit_mat(1,:),crit_mat2(1,:)] = compute_criteria(gd);

% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
    gabdualnorm(g,gd,a,M));
fprintf('  Equal to the canonical dual? norm(gd-gcan) = %g \n',...
    norm(gd-gcan));

figure();
plot_time(gd);
title('(c")');
save_name='l2_time';
plotfig(save_name,paramplot );
figure()
magresp(gd,'fir','dynrange',dr,'1');
title('(c")');
setplotfreq();
save_name='l2_freq';
plotfig(save_name,paramplot );



%% Convex optimization L1 norm in time

%parameter for convex optimization
alpha=0.01;% weight of the L1 norm in time
mu=0;      % smoothing parameter
gamma=0; % weight for the concentration in time
drawnow;
% solving the problem
fprintf(' Solve the convex optimization problem: L1 time\n');
gd=gaboptdual(g,a,M,'alpha',alpha,'mu',mu,'gamma',gamma,...
    'maxit',10*maxit,'tol',0.0001*tol,'evolution',evo,'constant',...
    'gif','l1_time');
gd=real(gd);

[crit_mat(6,:),crit_mat2(6,:)] = compute_criteria(gd);

% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
    gabdualnorm(g,gd,a,M));

figure();
plot_time(gd);
title('(a)');
save_name='l1time_time';
plotfig(save_name,paramplot );
figure()
magresp(gd,'fir','dynrange',dr,'1');
title('(b)');
setplotfreq();
save_name='l1time_freq';
plotfig(save_name,paramplot );

figure()
title('(a")');
ax=plot_time(gd);
hold off;
semilogy(ax,abs(fftshift(gd)));
ylim([1e-10 10*max(gd)]);
save_name='l1time_time_log';
plotfig(save_name,paramplot );
%%
% try to set to 0 coefficient below a certain threshold
thresh=80; %threshlod in dB
t=max(gd)*10^(-thresh/10);
gd(abs(gd)<t)=0;    %thresholding

% Compute reconstruciton error, number of non zeros coefficients and length
% of support
fprintf('  Reconstruction error of the window after a threshold of %i dB: %g \n',...
    thresh,gabdualnorm(g,gd,a,M));
fprintf('  Percent of non-zeros coefficients: %g \n',...
    sum(abs(gd)>0)/length(gd)*100);
fprintf('  Number of non-zeros coefficients: %g \n',...
    sum(abs(gd)>0));
fprintf('  Percent of non-zeros for the cannonical dual: %g \n',...
    sum(abs(gcan)>t)/length(gd)*100);

% length of the support
nz=fftshift(abs(gd))>0;
indi = find(nz, 1, 'first');
indf = find(nz, 1, 'last');
fprintf('  Length of the support: %i\n',indf-indi+1);

%% Convex optimization L1 norm in frequency

%parameter for convex optimization

beta=0.001; % weight of the L1 norm in frequency

mu=0;      % smoothing parameter
gamma=0; % weight for the concentration in time

% solving the problem
fprintf(' Solve the convex optimization problem: L1 frequency\n');
gd=gaboptdual(g,a,M,'beta',beta,'mu',mu,'gamma',gamma,...
    'maxit',maxit*10,'tol',0.01*tol,'evolution',evo,'constant',...
    'gif','l1_freq');
gd=real(gd);

[crit_mat(7,:),crit_mat2(7,:)] = compute_criteria(gd);

% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
    gabdualnorm(g,gd,a,M));

figure();
plot_time(gd);
title('(c)');
save_name='l1freq_time';
plotfig(save_name,paramplot );

figure()
magresp(gd,'fir','dynrange',dr,'1');
title('(d)');
setplotfreq();
save_name='l1freq_freq';
plotfig(save_name,paramplot );

%% Convex optimization Variance in time

%parameter for convex optimization
mu=0;      % smoothing parameter
gamma=0; % weight for the concentration in time
vart=1/100;
drawnow;
% solving the problem
fprintf(' Solve the convex optimization problem: Variance in time\n');
gd=gaboptdual(g,a,M,'mu',mu,'gamma',gamma,'vart',vart,...
    'maxit',3*maxit,'tol',0.1*tol,'evolution',evo,'dynamic',...
    'gif','var_time');
gd=real(gd);

[crit_mat(8,:),crit_mat2(8,:)] = compute_criteria(gd);

% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
    gabdualnorm(g,gd,a,M));

figure();
plot_time(gd);
title('(a)');
save_name='vart_time';
plotfig(save_name,paramplot );
figure()
magresp(gd,'fir','dynrange',dr,'1');
title('(b)');
setplotfreq();
save_name='vart_freq';
plotfig(save_name,paramplot );


%% Convex optimization: Energy variance in time

%parameter for convex optimization
mu=0;      % smoothing parameter
gamma=0; % weight for the concentration in time
var2t=10;
drawnow;
% solving the problem
fprintf(' Solve the convex optimization problem: Energy variance in time\n');
gd=gaboptdual(g,a,M,'mu',mu,'gamma',gamma,'var2t',var2t,...
    'maxit',maxit,'tol',tol,'evolution',evo,'constant',...
    'gif','var2_time');
gd=real(gd);

[crit_mat(10,:),crit_mat2(10,:)] = compute_criteria(gd);

% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
    gabdualnorm(g,gd,a,M));

figure();
plot_time(gd);
title('(c)');
save_name='var2t_time';
plotfig(save_name,paramplot );
figure()
magresp(gd,'fir','dynrange',dr,'1');
title('(d)');
setplotfreq();
save_name='var2t_freq';
plotfig(save_name,paramplot );

%% Convex optimization concentration in time

%parameter for convex optimization
mu=0;      % smoothing parameter
gamma=100; % weight for the concentration in time
drawnow;
% solving the problem
fprintf(' Solve the convex optimization problem: Concentration in time\n');
gd=gaboptdual(g,a,M,'mu',mu,...
    'gamma',gamma,'maxit',maxit,'tol',tol,'evolution',evo,'dynamic',...
    'gif','grad_freq');
gd=real(gd);

[crit_mat(3,:),crit_mat2(3,:)] = compute_criteria(gd);

% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
    gabdualnorm(g,gd,a,M));

figure();
plot_time(gd);
title('(e)');
save_name='smoothfreq_time';
plotfig(save_name,paramplot );
figure()
magresp(gd,'fir','dynrange',dr,'1');
title('(f)');
setplotfreq();
save_name='smoothfreq_freq';
plotfig(save_name,paramplot );

gd_time=gd;

%% Convex optimization Variance in frequency

%parameter for convex optimization
mu=0;      % smoothing parameter
gamma=0; % weight for the concentration in time
varf=8/100;
drawnow;

% solving the problem
fprintf(' Solve the convex optimization problem: Variance in frequency\n');
gd=gaboptdual(g,a,M,'mu',mu,'gamma',gamma,'varf',varf,...
   'maxit',10*maxit,'tol',0.01*tol,'evolution',evo,'dynamic',...
    'gif','var_freq');
gd=real(gd);

[crit_mat(9,:),crit_mat2(9,:)] = compute_criteria(gd);

% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
    gabdualnorm(g,gd,a,M));

figure();
plot_time(gd);
title('(a)');
save_name='varf_time';
plotfig(save_name,paramplot );
figure()
magresp(gd,'fir','dynrange',dr,'1');
title('(b)');
setplotfreq();
save_name='varf_freq';
plotfig(save_name,paramplot );


%% Convex optimization: Energy variance in frequency

%parameter for convex optimization
mu=0;      % smoothing parameter
gamma=0; % weight for the concentration in time
var2f=10;
drawnow;
% solving the problem
fprintf(' Solve the convex optimization problem: Energy variance in frequency\n');
gd=gaboptdual(g,a,M,'mu',mu,'gamma',gamma,'var2f',var2f,...
    'maxit',maxit,'tol',tol,'evolution',evo,'constant',...
    'gif','var2_freq');
gd=real(gd);

[crit_mat(11,:),crit_mat2(11,:)] = compute_criteria(gd);

% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
    gabdualnorm(g,gd,a,M));

figure();
plot_time(gd);
title('(c)');
save_name='var2f_time';
plotfig(save_name,paramplot );

figure()
magresp(gd,'fir','dynrange',dr,'1');
title('(d)');
setplotfreq();
save_name='var2f_freq';
plotfig(save_name,paramplot );

%% Convex optimization smoothness constraint

%parameter for convex optimization
gamma=0; % weight for the concentration in time
mu=1000;      % smoothing parameter

drawnow;
% solving the problem
fprintf(' Solve the convex optimization problem: smoothness constraint\n');
gd=gaboptdual(g,a,M,'gamma',gamma,'mu',mu,...
    'maxit',maxit,'tol',tol,'evolution',evo,'dynamic',...
    'gif','grad_time');
gd=real(gd);

[crit_mat(2,:),crit_mat2(2,:)] = compute_criteria(gd);

% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
    gabdualnorm(g,gd,a,M));
figure();
plot_time(gd);
title('(e)');
save_name='smoothtime_time';
plotfig(save_name,paramplot );
figure()
magresp(gd,'fir','dynrange',dr,'1');
title('(f)');
setplotfreq();
save_name='smoothtime_freq';
plotfig(save_name,paramplot );

gd_freq=gd;

%% Convex optimization S0 norm

%parameter for convex optimization
mu=0;      % smoothing parameter
gamma=0; % weight for the concentration in time
delta=0.005;

% solving the problem
fprintf(' Solve the convex optimization problem: S0 norm\n');
gd=gaboptdual(g,a,M,'mu',mu,'gamma',gamma,'delta',delta,...
    'maxit',maxit,'tol',tol,'quiet','evolution',evo,'dynamic',...
    'gif','S0');
gd=real(gd);

[crit_mat(4,:),crit_mat2(4,:)] = compute_criteria(gd);

% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
    gabdualnorm(g,gd,a,M));

figure();
plot_time(gd);
title('(a)');
save_name='S0_time';
plotfig(save_name,paramplot );
figure()
magresp(gd,'fir','dynrange',dr,'1');
title('(b)');
setplotfreq();
save_name='S0_freq';
plotfig(save_name,paramplot );

%% Convex optimization weighted S0 norm


%parameter for convex optimization
mu=0;      % smoothing parameter
gamma=0; % weight for the concentration in time
deltaw=0.005;

% solving the problem
fprintf(' Solve the convex optimization problem: S0 norm\n');
gd=gaboptdual(g,a,M,'mu',mu,'gamma',gamma,'deltaw',deltaw,...
    'maxit',maxit,'tol',tol,'quiet','evolution',evo,'dynamic',...
    'gif','S0_weighted');
gd=real(gd);

[crit_mat(5,:),crit_mat2(5,:)] = compute_criteria(gd);

% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
    gabdualnorm(g,gd,a,M));

figure();
plot_time(gd);
title('(c)');
save_name='weighted_S0_time';
plotfig(save_name,paramplot );
figure()
magresp(gd,'fir','dynrange',dr,'1');
title('(d)');
setplotfreq();
save_name='weighted_S0_freq';
plotfig(save_name,paramplot );



%% Convex optimization L2 - sine like

glike=sin(2*pi*[1:L]'/L);
% solving the problem
fprintf(' Solve the convex optimization problem: sine like\n');
gd=gabglike(g,glike,a,M,...
    'gif','sine_like','evolution',evo);
gd=real(gd);

% [crit_mat(12,:),crit_mat2(12,:)] = compute_criteria(gd);

% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
    gabdualnorm(g,gd,a,M));
figure();
plot_time(fir2long(glike,L)*norm(gcan)/norm(glike));
title('(a)');
save_name='sine_time';
plotfig(save_name,paramplot );
figure()
plot_time(gd);
title('(b)');
save_name='sinelike_time';
plotfig(save_name,paramplot );

%% Convex optimization L2 - dirac like

glike = zeros(size(g));
glike(1) = 1;
% solving the problem
fprintf(' Solve the convex optimization problem: dirac like\n');
gd=gabglike(g,glike,a,M,'maxit',1000,...
    'gif','dirac_like','evolution',evo);
gd=real(gd);

% [crit_mat(13,:),crit_mat2(13,:)] = compute_criteria(gd);

% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
    gabdualnorm(g,gd,a,M));
figure();
plot_time(fir2long(glike,L)*norm(gcan)/norm(glike));
title('(c)');
save_name='dirac_time';
plotfig(save_name,paramplot );
figure()
plot_time(gd);
title('(d)');
save_name='diraclike_time';
plotfig(save_name,paramplot );

figure()
plot_time(gd);
ylim([-0.5*max(gcan), 1.5*max(gcan)]);
title('(d`)');
save_name='diraclike_time_zoom';
plotfig(save_name,paramplot );


%%

Y = (crit_mat - repmat(min(crit_mat(1:11,:)),11,1)); Y = Y./repmat(max(Y(1:11,:)),11,1);
Y

crit_mat2