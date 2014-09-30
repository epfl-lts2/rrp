%RR_GABFIRTIGHT Optimization of dual gabor windows with support constraint
%                           
%   Reproducible research addendum for optimization of dual gabor windows
%   ---------------------------------------------------------------------
%   
%   DESIGNING GABOR WINDOWS USING CONVEX OPTIMIZATION
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
%   In order to use this matlab file you need the UNLocXbox toolbox. You
%   can download it on http://unlocbox.sourceforge.net and the LTFAT
%   toolbox. You can download it on http://ltfat.sourceforge.net
%   
%   The problem
%   -----------
%
%   As starting setup, we choose a Gabor system $G(g,30,60)$ with
%   an Itersine window of length $L_g = 60$. For this half-overlap,
%   redundancy $2$ situation, the Itersine window forms a tight, painless
%   frame with better joint TF concentration than other widely used
%   constructions for redundancy $2$ tight frames. We now attempt the
%   construction of a Gabor Parseval frame with redundancy $2$, using a
%   window function that further improves the TF concentration of the
%   Itersine window.  
%   
%   To gain some design freedom, we allow the tight window $h$ to have
%   support length $L_h = 360$. The gradient priors are used to promote a
%   window that is smooth and well-localized in both domains, leading
%   (formally) to the optimization problem       : 
%
%   .. h  = argmin_x   gamma || nabla F x||_2^2 + mu || nabla x ||_2^2
%
%   ..     such that    x is a tight windows and is compaclty supported
%
%   .. math:: \begin{split}  \text{gd}  =  \text{arg} \min_x   &  \gamma \| \nabla Fx\|_2^2 + \mu \| \nabla x \|_2^2 \\     \text{such that }& x \text{ is a tight windows }  \\     \text{and }& x \text{ is compaclty supported } \end{split}
%
%   with 
%         - $\mu$   :  Smoothing parameter
%         - $\gamma$:  Weight of the gradient in time
%         - $\mu$   :  Weight of the gradient in frequency
%         - $x$     :  Optimization variable
%         - $F$     :  The Fourier matrix 
%
%   As starting setup, we choose a Gabor system $G(g,30,60)$ with
%   an Itersine window of length $L_g = 60$. For this half-overlap,
%   redundancy $2$ situation, the Itersine window forms a tight, painless
%   frame with better joint TF concentration than other widely used
%   constructions for redundancy $2$ tight frames. We now attempt the
%   construction of a Gabor Parseval frame with redundancy $2$, using a
%   window function that further improves the TF concentration of the
%   Itersine window.  
%   
%   To gain some design freedom, we allow the tight window $g_t$ to have
%   support length $L_h = 360$. The gradient priors are used to promote a
%   window that is smooth and well-localized in both domains, leading
%   (formally) to the optimization problem.
%
%   For easier comparison, we tuned the result to have roughly the same
%   visual concentration in time. The result shown was obtained for the
%   regularization parameters $\gamma = 1$, $\mu = 5$ and shows
%   improved decay and side lobe attenuation, when compared to the
%   Itersine.
%
%   .. figure::
%
%      Results in time
%
%      
%
%   .. figure::
%
%      Results in frequency
%
%      
%
%   Since the problem is not convex, a good starting value and a good
%   timestep are crucial. We have obtained good results and dependable
%   convergence when choosing a starting window that is not too far from
%   what we aim for, i.e. it already has a good frame bound ratio $B/A$ for
%   the Gabor parameters $a,M$ and shows the properties we wish to promote
%   in the tight window, e.g. TF concentration. Note that, especially for
%   frames with small redundancy, it has been observed that a tradeoff
%   between localization and smoothness in TF exists between the analysis
%   and dual windows. Therefore, low redundancy Parseval frame windows
%   provide, in comparison, suboptimal TF concentration. 
%
%   In this particualr case, our starting point is the original itersine
%   window of length $60$.
%



%% initialization
clear;
close all;



global GLOBAL_save
global GLOBAL_baw
global GLOBAL_figpaper

% plotting
baw=GLOBAL_baw; % plot everything in black and white
dr = 100;
figpaper = GLOBAL_figpaper;
figpaper = 0;

if baw
    paramplot.pathfigure = 'figures/gabfirtight/'; 
else
    paramplot.pathfigure = 'figures_color/gabfirtight/'; 
end
paramplot.position = [100 100 300 200];
paramplot.titleweight ='bold';
paramplot.save=GLOBAL_save;
paramplot.baw=baw; % plot everything in black and white




%%
% parameters of frame and length
a=30;
M=60;
L=360;

% For testing the reconstruction
Llong=dgtlength(10*L,a,M);
Lplot=L;


% windows initial point
g=firwin('itersine',M);



%parameter for convex optimization
% those alpha should be smaller than the value of the windows !
mu=50;       % smoothing parameter
gamma=10;    % smoothing parameter in frequency
maxit=500; % maximum number of iteration
tol=1e-10;  % tolerance to stop iterating



gp=gabtight(g,a,M,L);


fprintf('  Reconstruction error of the original tight window: %g\n', ...
    gabdualnorm(gp,gp,a,M,Llong));


fprintf('Solve the optimization problem \n');


% call optimization routine
h=gabfirtight(L,rand(size(g)),a,M,'mu',mu,...
    'gamma',gamma,'maxit',maxit,'tol',tol,'quiet','evolution',1);


% Force the result to be real
h=real(h);

%%

if figpaper
    paramplot.position = [100 100 300 200];
    
    figure();
    plot_time(g);
    title('(a)');
    save_name='tight_time';
    plotfig(save_name,paramplot );
    figure();
    magresp(g,'fir','dynrange',dr,'1');
    title('(b)');
    setplotfreq();
    save_name='tight_freq';
    plotfig(save_name,paramplot );
   
    figure()
    plot_ambiguity(g,dr);
    title('(c)')
    save_name='tight_amb';
    plotfig(save_name,paramplot );
    
    figure();
    plot_time(h);
    title('(d)');
    save_name='opt_time';
    plotfig(save_name,paramplot );
    figure();
    magresp(h,'fir','dynrange',dr,'1');
    title('(e)');
    setplotfreq();
    save_name='opt_freq';
    plotfig(save_name,paramplot );
   
    figure()
    plot_ambiguity(h,dr);
    title('(f)')
    save_name='opt_amb';
    plotfig(save_name,paramplot );
    
else
    figure(1);
    plot_time(fir2long(gp,Lplot),'b');
    hold on
    plot_time(fir2long(h,Lplot),'r');
    if baw
        set_baw_color();
    end
    legend('Itersine','Optimization');
    title('(a)');
    save_name='time';
    paramplot.baw=0;
    plotfig(save_name,paramplot );
    paramplot.baw=baw;
    figure(2)
    magresp(gp,'fir','dynrange',dr,'1','opts',{'b'});
    hold on
    magresp(h,'fir','dynrange',dr,'1','opts',{'r'});
    if baw
        set_baw_color();
    end
    title('(b)');
    legend('Itersine','Optimization');

    setplotfreq();
    save_name='freq';
    paramplot.baw=0;
    plotfig(save_name,paramplot );
    paramplot.baw=baw;

end

% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
    gabdualnorm(h,h,a,M,Llong));

%%
Lcomp=dgtlength(2*L+1,a,M);
[A,B] = gabframebounds(g,a,M);

fac = 4;
                            
[crit_mat(1,:),crit_mat2(1,:)] = compute_criteria(g/sqrt(A),Lcomp,fac,L);
[crit_mat(2,:),crit_mat2(2,:)] = compute_criteria(h,Lcomp,fac,L);
crit_mat = crit_mat(:,[2,3,8,9,10,11,4,5,1,6,7]);


%%

XX = [100,100,1,1,10,10,1,1,10,1,1]'*[1,1];
crit = crit_mat.*XX';

XX2 = [100,100,1,1,1]'*[1,1];
crit2 = [crit_mat2.*XX2',crit(:,[1,2])];

%mat2tex(crit,'table')
mat2tex(crit2,'table')

