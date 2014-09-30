%RR_GABFIRDUAL Optimization of dual gabor windows with support constraint
%                           
%   Reproducible research addendum for optimization of dual gabor windows
%   ---------------------------------------------------------------------
%   
%   DESIGNING GABOR DUAL WINDOWS USING CONVEX OPTIMIZATION
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
%   The problem
%   -----------
%
%   In this this experiment, we will present a dual Gabor window with short
%   support in a nonpainless setup, i.e. for a system with few frequency
%   channels. This has previously been attempted by Stromer (see the paper
%   for more information). However, the solutions obtained by the
%   truncation method proposed by him are often badly localized in
%   frequency due to the tendency of the truncation method to yield
%   nonsmooth solutions, i.e.  solutions with 'jumps' or discontinuity-like
%   behavior in time. 
%
%   Here, we present the construction of a dual window that combines short
%   support with reasonable smoothness and therefore frequency
%   concentration. We start from a Gabor system $G(g,30,60)$ with
%   redundancy $2$. The analysis window $g$ is chosen as a Nuttall window
%   of length $L_g=120$ samples. 
%
%   .. figure::
%
%      Analysis window in time domain
%
%      The chosen windows is 'Nuttall'
%
%   .. figure::
%
%      Analysis window in frequency domain
%
%      The 'Nuttall' window is very well localized in the frequency domain.
%
%   We desire a dual window $h$ with the same support as $g$, i.e. $L_h =
%   L_g = 120$. Therefore, we set the following convex problem: 
%
%   .. h  = argmin_x   alpha ||x||_1 + beta ||F x||_1 
%   
%   ..               +  mu || nabla x ||_2^2 + gamma || nabla F x ||_2^2 
%
%   ..     such that    x is a dual windows of g and is compaclty supported
%
%   .. math:: \begin{split}  \text{h}  =  \text{arg} \min_x   & \alpha \|x\|_1 + \beta \|Fx\|_1 + \mu \| \nabla x \|_2^2 + \gamma \| \nabla \mathcal{F} x \|_2^2 \\     \text{such that }& x \text{ is a dual windows of }g  \\     \text{and }& x \text{ is compaclty supported } \end{split}
%
%   with 
%         - $\mu$   :  Smoothing parameter
%         - $\gamma$:  Concentration parameter
%         - $\alpha$:  Weight of the L1 norm in time
%         - $\beta$ :  Weight of the L1 norm in frequency
%         - $x$     :  Optimization variable
%         - $F$     :  The Fourier matrix 
%
%   We minimize the norm of the gradient of x in order to have a smooth
%   windows. The L1 norms help the dual windows to be concentrated in time
%   and frequency.
%
%   The results in Figure below (c,d) show the optimal dual window with
%   regards to the regularization parameters: $\alpha_1=\beta_{2}=0.001$
%   and $\mu_{3}=\gamma_{4}=1$. Those values were chosen experimentally.
%
%   .. figure::
%
%      Dual synthesis window in time domain
%
%      The optmain dual windows is quite smooth and is compactly supported.
%
%   .. figure::
%
%      Dual synthesis window in frequency domain
%
%      The dual windows is smooth. As a consequence it is confined in low
%      frequencies.
%
%   For comparison, we include the least-squares solution provided by the
%   truncation method below(e,f).
%
%   .. figure::
%
%      Dual synthesis window (Matrix pseudo inversion) in time domain
%
%      The optmain dual windows is not smooth but is compactly supported.
%
%   .. figure::
%
%      Dual synthesis window (Matrix pseudo inversion) in frequency domain
%
%      The dual windows is not really smooth. It is not really well
%      concentrated.
%
%   Minimizing the selected regularization functions improves upon the
%   desired features, in particular smoothness and localization by means of
%   minimizing the gradient. The addition of the sparsity terms is
%   supposed to suppress the tendency of the solution to have M-like
%   shape, i.e. multiple peaks. This is unwanted as it leads to windows
%   with ambiguous temporal or frequency position. Heuristically,
%   minimizing the $l_1$-norm pushes all big coefficients to similar
%   values, therefore achieving the suppression of multiple significant
%   peaks.
%
%   The solution provided is assumed to perform perfect reconstruction on
%   any signal with admissible length greater or equal to $L$. As in the
%   previous experiments, the maximum relative reconstruction error
%   ($4.5e^{-14}$) is of the order of the machine precision. 
%
%   To guarantee the canonical dual window to be supported on $L_h=L_g$,
%   we would be required to increase the number of frequency channels to
%   $M >= 120$, putting us in the painless case. Therefore the redundancy is
%   increased twofold, which is an unwanted side effect. Alternatively, in
%   this setting, we could decide to keep the parameters $a = 30$, $M = 60$
%   fixed, but decrease the window size to $L_g < =60$. However, this
%   construction provides a system with a more than $8$ times larger frame
%   bound ratio. Consequently, the resulting canonical dual window $h$,
%   shown in the figure bellow, shows bad frequency behavior and an
%   undesirable, M-like shape in time. In contrast, the method proposed in
%   the paper allows the use of nicely shaped,compactly supported dual
%   Gabor windows at low redundancies,  without the strong restrictions of
%   the painless case.   
%
%   .. figure::
%
%      Painless synthesis window in time domain
%
%      The painless window has also a bad time localisation
%
%   .. figure::
%
%      Painless synthesis window in frequency domain
%
%      The painless windows has a bad frequency localisation.
%





% initialization
clear;
close all;

global GLOBAL_save
global GLOBAL_baw



% parameters of frame and length
a=30;
M=60;
Lg=120;
Ldual=120;

% For testing the reconstruction
Llong=dgtlength(10*(Ldual+Lg),a,M);
Lplot=320;

% plotting
dr=100;
if GLOBAL_baw
    paramplot.pathfigure = 'figures/gabfirdual/'; 
else
    paramplot.pathfigure = 'figures_color/gabfirdual/'; 
end
paramplot.position = [100 100 300 200];
paramplot.titleweight ='bold';
paramplot.save = GLOBAL_save;
paramplot.baw = GLOBAL_baw; % plot everything in black and white

% windows
g=firwin('nuttall',Lg);
% g=randn(size(g));
% g=g/norm(g);


%parameter for convex optimization
% those alpha should be smaller than the value of the windows !
alpha=0.001; % weight of the L1 norm in time
beta=0.001; % weight of the L1 norm in frequency
mu=1;         % smoothing parameter in time
gamma=1;      % smoothing parameter in frequency

maxit=500;   % maximum number of iteration
tol=1e-6;    % tolerance to stop iterating




figure(1)
plot_time(fir2long(g,Lplot));
title('(a)');
save_name='nuttall_time';
plotfig(save_name,paramplot );

figure(2);
magresp(g,'fir','dynrange',dr,'1');
title('(b)');
setplotfreq();
save_name='nuttall_freq';
plotfig(save_name,paramplot );

fprintf('Solve the optimization problem \n');
% call optimization routine
h=gabfirdual(Ldual,g,a,M,'alpha',alpha,'beta',beta,'mu',mu,...
    'gamma',gamma,'maxit',maxit,'tol',tol,'quiet');

%gd_cvx=gabfirdual_cvx(Ldual,g,a,M,alpha,beta,mu,gamma);
%gd=gd_cvx;

h=real(h);

figure(3)
plot_time(fir2long(h,Lplot));
title('(c)');
save_name='sol_time';
plotfig(save_name,paramplot );

figure(4);
magresp(h,'fir','dynrange',dr,'1');
title('(d)');
setplotfreq();
save_name='sol_freq';
plotfig(save_name,paramplot );

% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
    gabdualnorm(g,h,a,M,Llong));


%%

[gd_pinv,Gcut]=gabfirdual_pinv(Ldual,g,a,M);
disp('Reconstruction error of the matrix inversion method.');
reconerr=gabdualnorm(g,gd_pinv,a,M,Llong);
fprintf('  Reconstruction error of the matrix pseudo-inversion method: %g \n',...
    gabdualnorm(g,gd_pinv,a,M,Llong));
gd_pinv=real(gd_pinv);

figure(5)
plot_time(fir2long(gd_pinv,Lplot));
title('(e)');
save_name='pinv_time';
plotfig(save_name,paramplot );

figure(6);
magresp(gd_pinv,'fir','dynrange',dr,'1');
title('(f)');
setplotfreq();
save_name='pinv_freq';
plotfig(save_name,paramplot );


%%

Lg=60;
g=firwin('nuttall',Lg);
gpainless=gabdual(g,a,M,Ldual);

figure(7)
plot_time(fir2long(gpainless,Lplot));
title('(a)');
save_name='painless_time';
plotfig(save_name,paramplot );

figure(8);
magresp(gpainless,'fir','dynrange',dr,'1');
title('(b)');
setplotfreq();
save_name='painless_freq';
plotfig(save_name,paramplot );


%%
Lcomp=dgtlength(Lg+Ldual+1,a,M);

fac = 4;
                            
[crit_mat(1,:),crit_mat2(1,:)] = compute_criteria(gd_pinv,Lcomp,fac,Ldual);
[crit_mat(2,:),crit_mat2(2,:)] = compute_criteria(h,Lcomp,fac,Ldual);
crit_mat = crit_mat(:,[2,3,8,9,10,11,4,5,1,6,7]);


%%
XX = [100,100,1,1,10,10,1,1,10,1,1]'*[1,1];
crit = crit_mat.*XX';

XX2 = [100,100,1,1,1]'*[1,1];
crit2 = [crit_mat2.*XX2',crit(:,[1,2,9,10])];

%mat2tex(crit,'table')
mat2tex(crit2,'table')
