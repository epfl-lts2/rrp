%RR_GABFIRDUAL_SMOOTHNESS3 Optimization of dual gabor windows with support constraint
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
%   This experiment is very similar to rr_gabfirdual_smoothness, but we
%   works with 10 times longer windows.
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
%   For this experiment, we have chosen an itersine window of length 
%   $L_g = 2400$ (see figure below , $a = 500$ and  $M=3000$. The support
%   of the dual window candidates was restricted to $L_h = 6000$.  
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
%   We apply 2 different methods for measuring joint time- frequency
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
%     We consider both the magnitude variance, as in Heisenberg's 
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
%   The figures below shows time and frequency representations, as well as
%   the ambiguity function of the results. We see that all the
%   criteria provide (visually) nicely concentrated dual Gabor windows,
%   improved over the canonical dual. In this case, we do not optimize the
%   S0 norm because it will be to computationally expensive.   
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
%   If you run the file, you'll notice that it takes a few second before
%   the optimization algorithm starts. This is due to the fact, that we
%   pseudo invert the WR matrix. Then each iteration is very quick.
%

%% Initialization
clear;
close all;
% clc


global GLOBAL_save
global GLOBAL_baw
global GLOBAL_figpaper
global GLOBAL_evolution

GLOBAL_evolution = 1;

if GLOBAL_evolution
    evo = GLOBAL_evolution;
else
    evo = 0;
end

% plotting
baw=GLOBAL_baw; % plot everything in black and white
dr = 160;
figpaper = GLOBAL_figpaper;
if baw
    paramplot.pathfigure = 'figures/gabfirdual_smooth_long/'; 
else
    paramplot.pathfigure = 'figures_color/gabfirdual_smooth_long/'; 
end
paramplot.position = [100 100 600 400];
paramplot.titleweight ='bold';
paramplot.save = GLOBAL_save;
paramplot.baw = baw; % plot everything in black and white


%%
% parameters of frame and length
a=500;
M=3000;
L=3000;


Ldual = 6000;
Ltest = 8*M;

% windows

g=fir2long(firwin('itersine',2400),L);

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

maxit=500;   % maximum number of iteration
tol=1e-4;    % tolerance to stop iterating




%% Solves the problems

%parameter for convex optimization
mu=10;         % smoothing parameter in time
gamma=10;      % smoothing parameter in frequency


fprintf('Solve the optimization problem : Smoothness\n');
% call optimization routine
gd=gabfirdual(Ldual,g,a,M,'maxit',maxit,'tol',tol,...
'mu',mu,'gamma',gamma,'evolution',evo,'debug');
gd_grad=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,Ltest));
%%
%parameter for convex optimization
mu=0;         % smoothing parameter in time
gamma=0;      % smoothing parameter in frequency
vart=3/L;   % parameter for the variance in time
varf=3/L;   % parameter for the variance in frequency


fprintf('Solve the optimization problem : Variance\n');
% call optimization routine
gd=gabfirdual(Ldual,g,a,M,'maxit', maxit,'tol',tol,...
'mu',mu,'gamma',gamma,'vart',vart,'varf',varf,'evolution',evo,'dynamic');
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
gd=gabfirdual(Ldual,g,a,M,'maxit',maxit,'tol',tol,...
'mu',mu,'gamma',gamma,'var2t',var2t,'var2f',var2f,'evolution',evo,'constant');
gd_var2=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,Ltest));


%%
if  figpaper
    
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
    plot_time(gcan);
    title('(c)');
    save_name='can_dual_time';
    plotfig(save_name,paramplot2 );
    figure(4);
    magresp(gcan,'fir','dynrange',dr,'1');
    title('(d)');
    setplotfreq();
    save_name='can_dual_freq';
    plotfig(save_name,paramplot2 );
   
    figure(5)
    plot_ambiguity(gcan,dr);
    title('(e)')
    save_name='can_dual_amb';
    plotfig(save_name,paramplot2 );

    figure(6);
    plot_time(gd_grad);
    title('(f)');
    save_name='grad_time';
    plotfig(save_name,paramplot2 );
    figure(7);
    magresp(gd_grad,'fir','dynrange',dr,'1');
    title('(g)');
    setplotfreq();
    save_name='grad_freq';
    plotfig(save_name,paramplot2 );
    
    figure(8)
    plot_ambiguity(gd_grad,dr);
    title('(h)')
    save_name='grad_amb';
    plotfig(save_name,paramplot2 );

    figure(9);
    plot_time(gd_var);
    title('(i)');
    save_name='var_time';
    plotfig(save_name,paramplot2 );
    figure(10);
    magresp(gd_var,'fir','dynrange',dr,'1');
    title('(j)');
    setplotfreq();
    save_name='var_freq';
    plotfig(save_name,paramplot2 );

    figure(11)
    plot_ambiguity(gd_var,dr);
    title('(k)')
    save_name='var_amb';
    plotfig(save_name,paramplot2 );
    
    figure(12);
    plot_time(gd_var2);
    title('(l)');
    save_name='var2_time';
    plotfig(save_name,paramplot2 );
    figure(13);
    magresp(gd_var2,'fir','dynrange',dr,'1');
    title('(m)');
    setplotfreq();
    save_name='var2_freq';
    plotfig(save_name,paramplot2 );

    figure(14)
    plot_ambiguity(gd_var2,dr);
    title('(n)')
    save_name='var2_amb';
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


    if baw
        set_baw_color();
    end

    legend('Canonical dual','Gradient','Variance','Energy Variance');
    title('Concentration in time-frequency');
    setplottime(gd_var);
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
    legend('Canonical dual','Gradient','Variance','Energy Variance');
    setplotfreq();
    save_name='comp_opt_freq';
    plotfig(save_name,paramplot );




    figure(5);
    paramplot.position = [100 100 600 600];
    save_name='comp_opt_ambiguity';
    dr = 160;
    subplot(221)
    plot_ambiguity(gcan,dr);
    title('Canonical dual');
    plotfig(save_name,paramplot );
    
    subplot(222)
    plot_ambiguity(gd_grad,dr);
    title('Gradient');
    plotfig(save_name,paramplot );
    
    subplot(223)
    plot_ambiguity(gd_var,dr);
    title('Variance');
    plotfig(save_name,paramplot );
    
    subplot(224)
    plot_ambiguity(gd_var2,dr);
    title('Energy variance');
    plotfig(save_name,paramplot );
    





    
end

%%
Lcomp=dgtlength(L+Ldual+1,a,M);

fac = 1;
[crit_mat(1,:),crit_mat2(1,:)] = compute_criteria(gcan,Lcomp,fac);
[crit_mat(2,:),crit_mat2(2,:)] = compute_criteria(gd_grad,Lcomp,fac);
[crit_mat(3,:),crit_mat2(3,:)] = compute_criteria(gd_var,Lcomp,fac);
[crit_mat(4,:),crit_mat2(4,:)] = compute_criteria(gd_var2,Lcomp,fac);

Y = (crit_mat(1:4,:) - repmat(min(crit_mat(1:4,:)),4,1)); 
Y = Y./repmat(max(Y(1:4,:)),4,1);

crit_mat = crit_mat(:,[2,3,8,9,10,11,4,5,1,6,7]);
%%
XX = [10000,1000,10,100,100,1000,1,1,100,10,10]'*[1,1,1,1];
crit = crit_mat.*XX'

XX2 = [100,100,1,1,1]'*[1,1,1,1];
crit2 = crit_mat2.*XX2'