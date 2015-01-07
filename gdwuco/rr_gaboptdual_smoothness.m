%RR_GABOPTDUAL_SMOOTHNESS Optimization of dual gabor windows with support constraint
%                           
%   Reproducible research addendum for optimization of dual gabor windows
%   ---------------------------------------------------------------------
%   
%   DESIGNING GABOR DUAL WINDOWS USING CONVEX OPTIMIZATION
%
%   Paper: Nathanael Perraudin, Nicki Holighaus, Peter L. Sondergaard and Peter Balazs
%   
%   Demonstration matlab file:  Perraudin Nathanael
%
%   ARI -- January 2013
%   
%   Dependencies
%   ------------
%
%   In order to use this matlab file you need the UNLocbox toolbox. You
%   can download it on https://lts2research.epfl.ch/unlocbox and the LTFAT
%   toolbox. You can download it on http://ltfat.sourceforge.net
%   
%   The problem
%   -----------
%
%   In this demo file, we solve the following problem: 
%
%   .. gd  = argmin_x    f(x)
%
%   ..     such that    x is a dual windows of g
%
%   .. math:: \begin{split}  \text{gd}  =  \text{arg} \min_x   & f(x) \\     \text{such that }& x \text{ is a dual windows of }g   \end{split}
%
%   with 
%         - $x$     :  Optimization variable
%         - $f(x)$  :  Functional to optimize
%
%   We present the construction of time-frequency concentrated dual Gabor 
%   windows for a system of fixed length $L = 500$. Our setup considers 
%   a frame $G(g,25,125)$, i.e. a system with redundancy $5$, where $g$ is 
%   a rectangular window of length $L_g = 125$ samples. This construction 
%   allows the comparison of alternative time-frequency localization 
%   measures and their effect on the solution window. The system forms a 
%   tight frame with badly localized dual window and is therefore an ideal 
%   sandbox for experimenting with various localization measures.
%
%   .. figure::
%
%      Analysis window in time domain
%
%      Rectangular tight window
%
%   .. figure::
%
%      Analysis windows in frequency domain
%
%      Rectangular tight window
%
%   We apply 3 different methods for measuring joint time- frequency
%   localization: 
%
%   * We combine the l2-norms of the gradient of $x$ and of $Fx$ (F beeing 
%     the Fourier transform), to optimize the smoothness in the time and the 
%     frequency domain. In that case:
%
%     ..  f(x) = ||nabla x||_2^2 + ||nabla Fx||_2^2
%
%     .. math:: f(x) = \|\nabla x\|_2^2 + \|\nabla Fx\|_2^2
%       
%   * Variance is associated with the Heisenberg uncertainty principle, a 
%     classical measure for time-frequency localiza- tion. While we cannot 
%     directly optimize the Heisenberg uncertainty (the product 
%     var(x) var(Fx) ), the sum of variances provides a very similar measure. 
%     We consider both the magnitude variance, as in Heisenberg?s 
%     uncertainty, as well as energy variance. In that case $f(x)$ can be
%     too different things
%
%     ..  f(x)= var(x) + var(Fx)
%
%     .. math:: f(x)= \text{var}(x) + \text{var}(Fx)
%
%     or
%
%     ..  f(x)= var(x^2) + var((Fx)^2)
%
%     .. math:: f(x)= \text{var}(x^2) + \text{var}((Fx)^2)
%     
%   * The S0-norm of $x$, or equivalently l1-norm of its STFT, is widely 
%     used in Gabor analysis to measure the spread of x in the 
%     time-frequency plane. Since truly sparse short-time Fourier 
%     transforms do not exist, it promotes functions that are concentrated 
%     around few points in the time-frequency plane. Since concentrated 
%     windows have small spread, they also have a small S0-norm. However, 
%     small spread is not a guarantee for concentration. To further improve 
%     the concentration, we also consider to weight the S0-norm as with the 
%     following:
%     
%     ..    W(f,t)= ln( 1 + x(t)^2 + y(f)^2 )
%
%     .. math::  W(f,t)= \sqrt{ x(t)^2 + y(f)^2 }
%     
%     with $x(i) = y(i)$ being the i-th entry of the vector
%
%     ..    1/sqrt(L) [-L/2,L/2+1,....,0,1,...,L/2].
%
%     .. math::  \frac{1}{\sqrt(L)} [-L/2,L/2+1,....,0,1,...,L/2].
%
%     In that case, we have:
%   
%     ..    f(x) = ||x||_S0 = || G x ||_1
%
%     .. math::  f(x) = \|x\|_{S_0} = \|Gx\|_1
%
%     with $G$ a full STFT with a Gaussian window. If the weighted case is
%     considered, we have:
%
%     ..    f(x) = || W G x ||_1
%
%     .. math::  f(x) =  \|WGx\|_1
%
%   Experiment 1
%   ------------
%
%   To show the effect of the gradient and variance measures when applied
%   in a single domain, this experiment shows the results obtained for 
%   optimizing the concentration in time only. Since the analysis system 
%   forms a tight frame, the canonical dual is a rectangular windows of 
%   length 125. In this example, magnitude variance minimization provides 
%   the strongest concentration in terms of support, albeit at the cost of 
%   the clear, slim peak obtained from the very similar energy variance and
%   gradient optimization.
%
%   .. figure::
%
%      Windows in time domain
%
%      Optimization of the concentration in time with different methods.
%
%   .. figure::
%
%      Windows in frequency domain
%
%      Optimization of the concentration in time with different methods.
%
%   .. figure::
%
%      Windows in time domain
%
%      Optimization of the concentration in frequency with different 
%      methods.
%
%   .. figure::
%
%      Windows in frequency domain
%
%      Optimization of the concentration in frequency with different 
%      methods.
%
%   Experiment 2
%   ------------
%
%   In fact, we are more interested in simultaneous time-frequency 
%   concentration. Therefore, we performed another experiment comparing the
%   time-frequency localization mea- sures proposed above. The following 
%   figures present the results of the experiment. The advantage of the 
%   gradient and variance approaches is the simplicity of tuning the 
%   tradeoff between time and frequency concentration, which involves only 
%   the appropriate choice of regularization parameter. On the other hand,
%   only S0-norm based optimization schemes truly measure concentration on
%   the whole time frequency plane, instead of its one dimensional 
%   projection on the time, respectively frequency domain.
%
%   The results show that pure S0-norm optimization is insufficient to 
%   obtain good concentration. Weighted S0-norm and magnitude variance 
%   optimization on the other hand, provide very good results. The 
%   magnitude variance shows slightly better concentration in time, while 
%   the weighted S0-norm solution is superior to all other measures in 
%   terms of time-frequency decay. As observed for the previous experiment,
%   energy variance and gradient approaches produce very similar results,
%   providing slightly worse time-frequency concentration compared to 
%   magnitude variance and weighted S0-norm. While the experiment suggests
%   weighted S0-norm and magnitude variance as the most promising 
%   approaches, this experiment should not be considered an exhaustive 
%   study of the various time-frequency localization measures. 
%
%   .. figure::
%
%      Windows in time domain
%
%      Optimization of the concentration in time-frequency with different 
%      methods.
%
%   .. figure::
%
%      Windows in frequency domain
%
%      Optimization of the concentration in time-frequency with different 
%      methods.
%
%   .. figure::
%
%      Spectrograms of the obtained windows
%
%      Optimization of the concentration in time-frequency with different 
%      methods.
%



%% initialization
clear all;
close all;





global GLOBAL_save
global GLOBAL_baw



% plotting
baw=GLOBAL_baw; % plot everything in black and white
dr=100;
if baw
    paramplot.pathfigure = 'figures/gaboptdual_smoothness/'; 
else
    paramplot.pathfigure = 'figures_color/gaboptdual_smoothness/'; 
end
paramplot.position = [100 100 600 400];
paramplot.titleweight ='bold';
paramplot.save = GLOBAL_save;
paramplot.baw = baw; % plot everything in black and white

%%

% parameters of frame and length
a=25;
M=125;
L=500;

% windows
g=fir2long(firwin('rect',L/4),L);
%g=pgauss(L,1);

% g=randn(size(g));
% g=g/norm(g);


%% Canonical dual
% The windows is tight. So the canonical dual is the same as the original
% window
fprintf('Compute the canonical dual \n');

gcan=gabdual(g,a,M,L);

%norm(gcan/norm(gcan)-g/norm(g))

% Verification: this number should be close to 0.
fprintf('  Reconstruction error of the canonical dual %g \n',...
                                gabdualnorm(g,gcan,a,M,L));


%%

maxit=1000;   % maximum number of iteration
tol=1e-5;    % tolerance to stop iterating


%%
paramplot2=paramplot;
paramplot2.position = [100 100 300 200];
figure(1);
plot_time(g);
title('(a)');
save_name='analysis_time';
plotfig(save_name,paramplot2 );
figure(2)
magresp(g,'fir','dynrange',dr,'1');
title('(b)');
setplotfreq();
save_name='analysis_freq';
plotfig(save_name,paramplot2 );


%% Solves the problems

%parameter for convex optimization
mu=0;         % smoothing parameter in time
gamma=1;      % smoothing parameter in frequency
delta=0.000;  % parameter for the S0 norm
deltaw=0.000; % parameter for the weighted norm  
vart=0/100;   % parameter for the variance in time
varf=0/100;   % parameter for the variance in frequency
var2t=0;      % parameter for the variance of the energy in time
var2f=0;      % parameter for the variance of the energy in frequency

fprintf('Solve the optimization problem : smoothness in frequency\n');
% call optimization routine
gd=gaboptdual(g,a,M,'maxit',maxit,'tol',tol,'delta',delta,'deltaw',deltaw,...
'mu',mu,'gamma',gamma,'vart',vart,'varf',varf,'var2t',var2t,'var2f',var2f);
gd_gradf=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,L));


%parameter for convex optimization
mu=0;         % smoothing parameter in time
gamma=0;      % smoothing parameter in frequency
delta=0.000;  % parameter for the S0 norm
deltaw=0.000; % parameter for the weighted norm  
vart=1/100;   % parameter for the variance in time
varf=0/100;   % parameter for the variance in frequency
var2t=0;      % parameter for the variance of the energy in time
var2f=0;      % parameter for the variance of the energy in frequency

fprintf('Solve the optimization problem : variance in time\n');
% call optimization routine
gd=gaboptdual(g,a,M,'maxit',maxit,'tol',tol,'delta',delta,'deltaw',deltaw,...
'mu',mu,'gamma',gamma,'vart',vart,'varf',varf,'var2t',var2t,'var2f',var2f);
gd_vart=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,L));

%parameter for convex optimization
mu=0;         % smoothing parameter in time
gamma=0;      % smoothing parameter in frequency
delta=0.000;  % parameter for the S0 norm
deltaw=0.000; % parameter for the weighted norm  
vart=0/100;   % parameter for the variance in time
varf=0/100;   % parameter for the variance in frequency
var2t=1;      % parameter for the variance of the energy in time
var2f=0;      % parameter for the variance of the energy in frequency

fprintf('Solve the optimization problem : energy variance in time\n');
% call optimization routine
gd=gaboptdual(g,a,M,'maxit',maxit,'tol',tol,'delta',delta,'deltaw',deltaw,...
'mu',mu,'gamma',gamma,'vart',vart,'varf',varf,'var2t',var2t,'var2f',var2f);
gd_var2t=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,L));

%% Display results
figure(3);
hold on;
plot_time(gcan,'k');
plot_time(gd_gradf,'b');
plot_time(gd_vart,'r');
plot_time(gd_var2t,'g');
if baw
    set_baw_color();
end
legend('Canonical dual','Gradient in freq','Variance','Energy Variance');
title('Concentration in time');
setplottime(gd);
xlim([0,80]);
save_name='comp_opt_time_time';
paramplot.baw=0;
plotfig(save_name,paramplot );
paramplot.baw=baw;
figure(4)
hold on;
magresp(gcan,'fir','dynrange',dr,'1','opts',{'k'});
magresp(gd_gradf,'fir','dynrange',dr,'1','opts',{'b'});
magresp(gd_vart,'fir','dynrange',dr,'1','opts',{'r'});
magresp(gd_var2t,'fir','dynrange',dr,'1','opts',{'g'});
legend('Canonical dual','Gradient in freq','Variance','Energy Variance');
title('(d)');
setplotfreq();
save_name='comp_opt_time_freq';
plotfig(save_name,paramplot );



%% Solves the problems

%parameter for convex optimization
mu=1;         % smoothing parameter in time
gamma=0;      % smoothing parameter in frequency
delta=0.000;  % parameter for the S0 norm
deltaw=0.000; % parameter for the weighted norm  
vart=0/100;   % parameter for the variance in time
varf=0/100;   % parameter for the variance in frequency
var2t=0;      % parameter for the variance of the energy in time
var2f=0;      % parameter for the variance of the energy in frequency

fprintf('Solve the optimization problem : smoothness in time\n');
% call optimization routine
gd=gaboptdual(g,a,M,'maxit',maxit,'tol',tol,'delta',delta,'deltaw',deltaw,...
'mu',mu,'gamma',gamma,'vart',vart,'varf',varf,'var2t',var2t,'var2f',var2f);
gd_gradt=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,L));


%parameter for convex optimization
mu=0;         % smoothing parameter in time
gamma=0;      % smoothing parameter in frequency
delta=0.000;  % parameter for the S0 norm
deltaw=0.000; % parameter for the weighted norm  
vart=0/100;   % parameter for the variance in time
varf=1/100;   % parameter for the variance in frequency
var2t=0;      % parameter for the variance of the energy in time
var2f=0;      % parameter for the variance of the energy in frequency

fprintf('Solve the optimization problem : variance in frequency\n');
% call optimization routine
gd=gaboptdual(g,a,M,'maxit',maxit,'tol',tol,'delta',delta,'deltaw',deltaw,...
'mu',mu,'gamma',gamma,'vart',vart,'varf',varf,'var2t',var2t,'var2f',var2f);
gd_varf=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,L));

%parameter for convex optimization
mu=0;         % smoothing parameter in time
gamma=0;      % smoothing parameter in frequency
delta=0.000;  % parameter for the S0 norm
deltaw=0.000; % parameter for the weighted norm  
vart=0/100;   % parameter for the variance in time
varf=0/100;   % parameter for the variance in frequency
var2t=0;      % parameter for the variance of the energy in time
var2f=1;      % parameter for the variance of the energy in frequency

fprintf('Solve the optimization problem : energy variance in frequency\n');
% call optimization routine
gd=gaboptdual(g,a,M,'maxit',maxit,'tol',tol,'delta',delta,'deltaw',deltaw,...
'mu',mu,'gamma',gamma,'vart',vart,'varf',varf,'var2t',var2t,'var2f',var2f);
gd_var2f=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,L));

%%
figure(5);
hold on;
plot_time(gcan,'k');
plot_time(gd_gradt,'b');
plot_time(gd_varf,'r');
plot_time(gd_var2f,'g');
legend('Canonical dual','Gradient in time','Variance in frequency',...
    'Energy Variance in frequency');
title('(e)');
setplottime(gd);
save_name='comp_opt_freq_time';
plotfig(save_name,paramplot );
figure(6)
hold on;
magresp(gcan,'fir','dynrange',dr,'1','opts',{'k'});
magresp(gd_gradt,'fir','dynrange',dr,'1','opts',{'b'});
magresp(gd_varf,'fir','dynrange',dr,'1','opts',{'r'});
magresp(gd_var2f,'fir','dynrange',dr,'1','opts',{'g'});
if baw
    set_baw_color();
end
legend('Canonical dual','Gradient in time','Variance in frequency',...
    'Energy Variance in frequency');
title('Concentration in frequency');
setplotfreq();
save_name='comp_opt_freq_freq';
paramplot.baw=0;
plotfig(save_name,paramplot );
paramplot.baw=baw;



%% Solves the problems

%parameter for convex optimization
mu=1;         % smoothing parameter in time
gamma=1;      % smoothing parameter in frequency
delta=0.000;  % parameter for the S0 norm
deltaw=0.000; % parameter for the weighted norm  
vart=0/100;   % parameter for the variance in time
varf=0/100;   % parameter for the variance in frequency
var2t=0;      % parameter for the variance of the energy in time
var2f=0;      % parameter for the variance of the energy in frequency

fprintf('Solve the optimization problem : Smoothness\n');
% call optimization routine
gd=gaboptdual(g,a,M,'maxit',maxit,'tol',tol,'delta',delta,'deltaw',deltaw,...
'mu',mu,'gamma',gamma,'vart',vart,'varf',varf,'var2t',var2t,'var2f',var2f);
gd_grad=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,L));

%parameter for convex optimization
mu=0;         % smoothing parameter in time
gamma=0;      % smoothing parameter in frequency
delta=0.000;  % parameter for the S0 norm
deltaw=0.000; % parameter for the weighted norm  
vart=1/100;   % parameter for the variance in time
varf=1/100;   % parameter for the variance in frequency
var2t=0;      % parameter for the variance of the energy in time
var2f=0;      % parameter for the variance of the energy in frequency

fprintf('Solve the optimization problem : Variance\n');
% call optimization routine
gd=gaboptdual(g,a,M,'maxit',maxit,'tol',tol,'delta',delta,'deltaw',deltaw,...
'mu',mu,'gamma',gamma,'vart',vart,'varf',varf,'var2t',var2t,'var2f',var2f);
gd_var=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,L));
                            
                            

%parameter for convex optimization
mu=0;         % smoothing parameter in time
gamma=0;      % smoothing parameter in frequency
delta=0.000;  % parameter for the S0 norm
deltaw=0.000; % parameter for the weighted norm  
vart=0/100;   % parameter for the variance in time
varf=0/100;   % parameter for the variance in frequency
var2t=1;      % parameter for the variance of the energy in time
var2f=1;      % parameter for the variance of the energy in frequency

fprintf('Solve the optimization problem : energy variance\n');
% call optimization routine
gd=gaboptdual(g,a,M,'maxit',maxit,'tol',tol,'delta',delta,'deltaw',deltaw,...
'mu',mu,'gamma',gamma,'vart',vart,'varf',varf,'var2t',var2t,'var2f',var2f);
gd_var2=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,L));


%parameter for convex optimization
mu=0;         % smoothing parameter in time
gamma=0;      % smoothing parameter in frequency
delta=0.001;  % parameter for the S0 norm
deltaw=0.000; % parameter for the weighted norm  
vart=0/100;   % parameter for the variance in time
varf=0/100;   % parameter for the variance in frequency
var2t=0;      % parameter for the variance of the energy in time
var2f=0;      % parameter for the variance of the energy in frequency

fprintf('Solve the optimization problem : S0\n');
% call optimization routine
gd=gaboptdual(g,a,M,'maxit',maxit,'tol',tol,'delta',delta,'deltaw',deltaw,...
'mu',mu,'gamma',gamma,'vart',vart,'varf',varf,'var2t',var2t,'var2f',var2f);
gd_S0=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,L));

%parameter for convex optimization
mu=0;         % smoothing parameter in time
gamma=0;      % smoothing parameter in frequency
delta=0.000;  % parameter for the S0 norm
deltaw=0.001; % parameter for the weighted norm  
vart=0/100;   % parameter for the variance in time
varf=0/100;   % parameter for the variance in frequency
var2t=0;      % parameter for the variance of the energy in time
var2f=0;      % parameter for the variance of the energy in frequency

fprintf('Solve the optimization problem : Weigted S0\n');
% call optimization routine
gd=gaboptdual(g,a,M,'maxit',maxit,'tol',tol,'delta',delta,'deltaw',deltaw,...
'mu',mu,'gamma',gamma,'vart',vart,'varf',varf,'var2t',var2t,'var2f',var2f);
gd_S0w=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,L));
                            
%%                          
%parameter for convex optimization
mu=0;         % smoothing parameter in time
gamma=0;      % smoothing parameter in frequency
deltaw2=1; % parameter for the weighted norm  


fprintf('Solve the optimization problem : Weigted S2 norm\n');
% call optimization routine
gd=gaboptdual(g,a,M,'maxit',maxit,'tol',tol,'deltaw2',deltaw2,...
'mu',mu,'gamma',gamma,'evolution',1,'debug');
gd_S0w2=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,L));

%%
figure(7);
hold on;
plot_time(gcan,'c');
plot_time(gd_grad,'b');
plot_time(gd_var,'r');
plot_time(gd_var2,'g');
plot_time(gd_S0,'k');
plot_time(gd_S0w,'m');
plot_time(gd_S0w2,'y');

if baw
    set_baw_color();
end

legend('Canonical dual','Gradient','Variance','Energy Variance',...
    'S0 norm','Weighted S0 norm','Weighted S2 norm');
title('Concentration in time-frequency');
setplottime(gd_var);
xlim([0,80]);
save_name='comp_opt_time';
paramplot.baw=0;
plotfig(save_name,paramplot );
paramplot.baw=baw;
figure(8)
hold on;
magresp(gcan,'fir','dynrange',dr,'1','opts',{'c'});
magresp(gd_grad,'fir','dynrange',dr,'1','opts',{'b'});
magresp(gd_var,'fir','dynrange',dr,'1','opts',{'r'});
magresp(gd_var2,'fir','dynrange',dr,'1','opts',{'g'});
magresp(gd_S0,'fir','dynrange',dr,'1','opts',{'k'});
magresp(gd_S0w,'fir','dynrange',dr,'1','opts',{'m'});
magresp(gd_S0w2,'fir','dynrange',dr,'1','opts',{'y'});
legend('Canonical dual','Gradient','Variance','Energy Variance',...
    'S0 norm','Weighted S0 norm','Weighted S2 norm');title('(h)');
setplotfreq();
save_name='comp_opt_freq';
plotfig(save_name,paramplot );

%%
figure(9);
paramplot.position = [100 100 600 600];
save_name='comp_opt_sgram';

subplot(321)
sgram(g/max(abs(g)),1,L,'dynrange',100);
title('Canonical dual');
plotfig(save_name,paramplot );

subplot(322)
sgram(gd_grad/max(abs(gd_grad)),1,L,'dynrange',100);
title('Gradient');
plotfig(save_name,paramplot );

subplot(323)
sgram(gd_var/max(abs(gd_var)),1,L,'dynrange',100);
title('Variance');
plotfig(save_name,paramplot );

subplot(324)
sgram(gd_var2/max(abs(gd_var2)),1,L,'dynrange',100);
title('Energy variance');
plotfig(save_name,paramplot );

subplot(325)
sgram(gd_S0/max(abs(gd_S0)),1,L,'dynrange',100);
title('S0 norm')
plotfig(save_name,paramplot );

subplot(326)
sgram(gd_S0w/max(abs(gd_S0w)),1,L,'dynrange',100);
title('Weighted S0 norm')
plotfig(save_name,paramplot );
                            
figure
sgram(gd_S0w2/max(abs(gd_S0w2)),1,L,'dynrange',100);
title('Weighted S2 norm')     

%%
crit_mat(1,:) = compute_criteria(gd);
crit_mat(2,:) = compute_criteria(gd_grad);
crit_mat(3,:) = compute_criteria(gd_var);
crit_mat(4,:) = compute_criteria(gd_var2);
crit_mat(5,:) = compute_criteria(gd_S0);
crit_mat(6,:) = compute_criteria(gd_S0w);
Y = (crit_mat - repmat(min(crit_mat(1:6,:)),6,1)); Y = Y./repmat(max(Y(1:6,:)),6,1);
Y
