%RR_GABFIRDUAL_SMOOTHNESS Optimization of dual gabor windows with support constraint
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
%   can download it on http://unlocbox.sourceforge.net and the LTFAT
%   toolbox. You can download it on http://ltfat.sourceforge.net
%   
%   The problem
%   -----------
%
%   In this experiment, we simultaneously optimize TF concentration of the
%   dual window. This is either achieved by a single prior on the TF
%   representation of the dual window $h$ (S0 and weighted S0 norm),
%   or by applying the time and frequency  versions of one prior, with
%   equal weights (gradient and variance). Generally, we write our problem:
%
%   .. gd  = argmin_x    f(x)
%
%   ..     such that    x is a dual windows of g and compactly supported
%
%   .. math:: \begin{split}  \text{gd}  =  \text{arg} \min_x   & f(x) \\     \text{such that }& x \text{ is a dual windows of }g \\     \text{such that }& x \text{ is compactly supported}  \end{split}
%
%   with 
%         - $x$     :  Optimization variable
%         - $f(x)$  :  Functional to optimize.
% 
%   For this experiment, we have chosen a Tukey window with a transition area
%   ratio of $3/5$ (see figure below, with $L_g = 240$, $a = 50$ and
%   $M=300$. The support of the dual window candidates was restricted to
%   $L_h = 600$. 
%
%   .. figure::
%
%      Analysis window in time domain
%
%      
%
%   .. figure::
%
%      Analysis windows in frequency domain
%
%      
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
%   The figures below shows time and frequency representations, as well as
%   the ambiguity function of the results. We see that all the
%   criteria provide (visually) nicely concentrated dual Gabor windows,
%   improved over the canonical dual. In particular, we see that the
%   gradient and energy variance optimal windows are very similar around
%   the origin, whereas the latter shows worse decay. The variance induces
%   the best concentration of large values in an almost rectangular TF
%   area, but similar to the energy variance, no good decay is achieved.
%   Both $S_0$ and weighted $S_0$ priors perform very well, but the weight
%   induces a more symmetric TF concentration and slightly better decay.
%   
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
%      Ambiguity functions of the obtained windows
%
%      Optimization of the concentration in time-frequency with different 
%      methods.
%
%   Note that, while the shape of this experiment's results is quite
%   characteristic, their quality cannot be representative for any
%   arbitrary setup. In fact the results are highly dependent on the
%   quality of the original Gabor frame, its redundancy and the choice of
%   support constraint.
%


%% Initialization
clear;
close all;



global GLOBAL_save
global GLOBAL_baw
global GLOBAL_figpaper

% plotting
baw=GLOBAL_baw; % plot everything in black and white
dr = 160;
figpaper = GLOBAL_figpaper;

if baw
    paramplot.pathfigure = 'figures/gabfirdual_smooth_tukeywin_06/'; 
else
    paramplot.pathfigure = 'figures_color/gabfirdual_smooth_tukeywin_06/'; 
end
paramplot.position = [100 100 600 400];
paramplot.titleweight ='bold';
paramplot.save = GLOBAL_save;
paramplot.baw = baw; % plot everything in black and white


%%
% parameters of frame and length
a=50;
M=300;
L=300;


Ldual = 600;
Ltest = 8*M;

% windows
%g=fir2long(firwin('rect',L/4),L);
 g = fir2long(fftshift(tukeywin(240,.6)),L);
% g=fir2long(firwin('hann',240),L);
% g=fir2long(firwin('itersine',240),L);

%g=pgauss(L,1);

% g=randn(size(g));
% g=g/norm(g);


%% Canonical dual
% The windows is tight. So the canonical dual is the same as the original
% window
fprintf('Compute the canonical dual \n');
if Ldual>M
gcan=gabdual(g,a,M,Ldual);
else
gcan=long2fir(gabdual(g,a,M,L),Ldual);
end



% Verification: this number should be close to 0.
fprintf('  Reconstruction error of the canonical dual %g \n',...
                                gabdualnorm(g,gcan,a,M,Ltest));


%%

maxit=150;   % maximum number of iteration
tol=1e-3;    % tolerance to stop iterating




%% Solves the problems

%parameter for convex optimization
mu=100;         % smoothing parameter in time
gamma=100;      % smoothing parameter in frequency


fprintf('Solve the optimization problem : Smoothness\n');
% call optimization routine
gd=gabfirdual(Ldual,g,a,M,'maxit',4*maxit,'tol',0.01*tol,...
'mu',mu,'gamma',gamma,'evolution',0);
gd_grad=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,Ltest));
%%
%parameter for convex optimization
mu=0;         % smoothing parameter in time
gamma=0;      % smoothing parameter in frequency
vart=30/L;   % parameter for the variance in time
varf=30/L;   % parameter for the variance in frequency


fprintf('Solve the optimization problem : Variance\n');
% call optimization routine
gd=gabfirdual(Ldual,g,a,M,'maxit', 4*maxit,'tol',0.01*tol,...
'mu',mu,'gamma',gamma,'vart',vart,'varf',varf,'evolution',0,'dynamic');
gd_var=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,Ltest));
 %%                           
                            

%parameter for convex optimization
mu=0;         % smoothing parameter in time
gamma=0;      % smoothing parameter in frequency
var2t=5;      % parameter for the variance of the energy in time
var2f=5;      % parameter for the variance of the energy in frequency

fprintf('Solve the optimization problem : energy variance\n');
% call optimization routine
gd=gabfirdual(Ldual,g,a,M,'maxit',4*maxit,'tol',0.01*tol,...
'mu',mu,'gamma',gamma,'var2t',var2t,'var2f',var2f,'evolution',0,'dynamic');
gd_var2=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,Ltest));
%%

%parameter for convex optimization
mu=0;         % smoothing parameter in time
gamma=0;      % smoothing parameter in frequency
delta=0.0003;  % parameter for the S0 norm


fprintf('Solve the optimization problem : S0\n');
% call optimization routine
gd=gabfirdual(Ldual,g,a,M,'maxit',maxit,'tol',tol,'delta',delta,...
'mu',mu,'gamma',gamma,'evolution',0,'debug','constant');
gd_S0=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,Ltest));
%%
%parameter for convex optimization
mu=0;         % smoothing parameter in time
gamma=0;      % smoothing parameter in frequency
deltaw=0.001; % parameter for the weighted norm  

fprintf('Solve the optimization problem : Weigted S0\n');
% call optimization routine
gd=gabfirdual(Ldual,g,a,M,'maxit',10*maxit,'tol',0.1*tol,'deltaw',deltaw,...
'mu',mu,'gamma',gamma,'evolution',0,'debug','dynamic');
gd_S0w=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,Ltest));

%%
if figpaper

    paramplot2=paramplot;
    paramplot2.position = [100 100 220 165];

    g = fir2long(g,Ldual);
    gcan = fir2long(gcan,Ldual);
    figure(1);
    plot_time(g);
    title('(a)');
    save_name='analysis_time';
    
    plotfig(save_name,paramplot2 );
    figure(2);
    magresp(g,'fir','dynrange',dr,'1');
    title('(b)');
    setplotfreq();
    save_name='analysis_freq';
    plotfig(save_name,paramplot2 );

    figure(3);
    plot_time(gcan);
    title('(a)');
    ylatex('||x||_2^2');
    
    save_name='can_dual_time';
    plotfig(save_name,paramplot2 );
  
    figure(4);
    magresp(gcan,'fir','dynrange',dr,'1');
    title('(b)');
    setplotfreq();
    save_name='can_dual_freq';
    plotfig(save_name,paramplot2 );
   
    figure(5)
    plot_ambiguity(gcan,dr);
    title('(c)')
    save_name='can_dual_amb';
    plotfig(save_name,paramplot2 );

    figure(6);
    plot_time(gd_grad);
    title('(d)');
    ylatex('||\nabla x||_2^2 + ||\nabla F x||_2^2');
    save_name='grad_time';
    plotfig(save_name,paramplot2 );
    figure(7);
    magresp(gd_grad,'fir','dynrange',dr,'1');
    title('(e)');
    setplotfreq();
    save_name='grad_freq';
    plotfig(save_name,paramplot2 );
    
    figure(8)
    plot_ambiguity(gd_grad,dr);
    title('(f)')
    save_name='grad_amb';
    plotfig(save_name,paramplot2 );

    figure(9);
    plot_time(gd_var);
    ylatex('var(x) + var(F x)');
    title('(g)');
    save_name='var_time';
    plotfig(save_name,paramplot2 );
    figure(10);
    magresp(gd_var,'fir','dynrange',dr,'1');
    title('(h)');
    setplotfreq();
    save_name='var_freq';
    plotfig(save_name,paramplot2 );

    figure(11)
    plot_ambiguity(gd_var,dr);
    title('(i)')
    save_name='var_amb';
    plotfig(save_name,paramplot2 );
    
    figure(12);
    plot_time(gd_var2);
    title('(j)'); 
    ylatex('var(x^2) + var((F x)^2)');
    save_name='var2_time';
    plotfig(save_name,paramplot2 );
    figure(13);
    magresp(gd_var2,'fir','dynrange',dr,'1');
    title('(k)');
    setplotfreq();
    save_name='var2_freq';
    plotfig(save_name,paramplot2 );

    figure(14)
    plot_ambiguity(gd_var2,dr);
    title('(l)')
    save_name='var2_amb';
    plotfig(save_name,paramplot2 );
    
    figure(15);
    plot_time(gd_S0);
    title('(m)');
    ylatex('||x||_{S0}');
    save_name='S0_time';
    plotfig(save_name,paramplot2 );
    figure(16);
    magresp(gd_S0,'fir','dynrange',dr,'1');
    title('(n)');
    setplotfreq();
    save_name='S0_freq';
    plotfig(save_name,paramplot2 );

    figure(17)
    plot_ambiguity(gd_S0,dr);
    title('(o)')
    save_name='S0_amb';
    plotfig(save_name,paramplot2 );
    
    figure(18);
    plot_time(gd_S0w);
    title('(p)');
    ylatex('||x||_{S0, weighted}');
    save_name='weighted_S0_time';
    plotfig(save_name,paramplot2 );
    figure(19);
    magresp(gd_S0w,'fir','dynrange',dr,'1');
    title('(q)');
    setplotfreq();
    save_name='weighted_S0_freq';
    plotfig(save_name,paramplot2 );

    figure(20)
    plot_ambiguity(gd_S0w,dr);
    title('(r)')
    save_name='weighted_S0_amb';
    plotfig(save_name,paramplot2 );
    
else
    
    paramplot2=paramplot;
    paramplot2.position = [100 100 300 200];
    
    g = fir2long(g,Ldual);
    gcan = fir2long(gcan,Ldual);
    figure(1);
    plot_time(g);
    title('(a)');
    save_name='analysis_time';
    plotfig(save_name,paramplot2 );
    figure(2);
    magresp(g,'fir','dynrange',dr,'1');
    title('(b)');
    setplotfreq();
    save_name='analysis_freq';
    plotfig(save_name,paramplot2 );
    
    figure(3);
    hold on;
    plot_time(gcan,'c');
    plot_time(gd_grad,'b');
    plot_time(gd_var,'r');
    plot_time(gd_var2,'g');
    plot_time(gd_S0,'k');
    plot_time(gd_S0w,'m');

    if baw
        set_baw_color();
    end

    legend('Canonical dual','Gradient','Variance','Energy Variance',...
        'S0 norm','Weighted S0 norm');
    title('Concentration in time-frequency');
    setplottime(gd_var);
    xlim([0,300]);
    save_name='comp_opt_time';
    paramplot.baw=0;
    plotfig(save_name,paramplot );
    paramplot.baw=baw;
    figure(4)
    hold on;
    magresp(gcan,'fir','dynrange',dr,'1','opts',{'c'});
    magresp(gd_grad,'fir','dynrange',dr,'1','opts',{'b'});
    magresp(gd_var,'fir','dynrange',dr,'1','opts',{'r'});
    magresp(gd_var2,'fir','dynrange',dr,'1','opts',{'g'});
    magresp(gd_S0,'fir','dynrange',dr,'1','opts',{'k'});
    magresp(gd_S0w,'fir','dynrange',dr,'1','opts',{'m'});
    legend('Canonical dual','Gradient','Variance','Energy Variance',...
        'S0 norm','Weighted S0 norm');title('(h)');
    setplotfreq();
    save_name='comp_opt_freq';
    plotfig(save_name,paramplot );




    figure(5);
    paramplot.position = [100 100 600 600];
    save_name='comp_opt_ambiguity';
    dr = 160;
    subplot(321)
    plot_ambiguity(gcan,dr);
    title('Canonical dual');
    plotfig(save_name,paramplot );
    
    subplot(322)
    plot_ambiguity(gd_grad,dr);
    title('Gradient');
    plotfig(save_name,paramplot );
    
    subplot(323)
    plot_ambiguity(gd_var,dr);
    title('Variance');
    plotfig(save_name,paramplot );
    
    subplot(324)
    plot_ambiguity(gd_var2,dr);
    title('Energy variance');
    plotfig(save_name,paramplot );
    
    subplot(325)
    plot_ambiguity(gd_S0,dr);
    title('S0 norm')
    plotfig(save_name,paramplot );
    
    subplot(326)
    plot_ambiguity(gd_S0w,dr);
    title('Weighted S0 norm')
    plotfig(save_name,paramplot );

end

%%
Lcomp=dgtlength(L+Ldual+1,a,M);

fac = 4;


[crit_mat(1,:),crit_mat2(1,:)] = compute_criteria(gcan,Lcomp,fac,Ldual);
[crit_mat(2,:),crit_mat2(2,:)] = compute_criteria(gd_grad,Lcomp,fac,Ldual);
[crit_mat(3,:),crit_mat2(3,:)] = compute_criteria(gd_var,Lcomp,fac,Ldual);
[crit_mat(4,:),crit_mat2(4,:)] = compute_criteria(gd_var2,Lcomp,fac,Ldual);
[crit_mat(5,:),crit_mat2(5,:)] = compute_criteria(gd_S0,Lcomp,fac,Ldual);
[crit_mat(6,:),crit_mat2(6,:)] = compute_criteria(gd_S0w,Lcomp,fac,Ldual);
%%
Y = (crit_mat(1:6,:) - repmat(min(crit_mat(1:6,:)),6,1)); 
Y = Y./repmat(max(Y(1:6,:)),6,1);
crit_mat = crit_mat(:,[2,3,8,9,10,11,4,5,1,6,7]);
%%
XX = [1000,1000,100,100,1000,1000,1,1,100,10,10]'*[1,1,1,1,1,1];
crit = crit_mat.*XX'

XX2 = [100,100,1,1,1]'*[1,1,1,1,1,1];
crit2 = crit_mat2.*XX2'

mat2tex(crit,'table');
mat2tex(crit2,'table');