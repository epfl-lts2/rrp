%RR_BEAT_ITERSINE Optimization of dual gabor windows with support constraint
%                           
%   Reproducible research addendum for optimization of dual gabor windows
%   ---------------------------------------------------------------------
%   
%   DESINGING GABOR DUAL WINDOWS USING CONVEX OPTIMIZATION
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
%   For low redudancy setting, there is usually a tradeoff between
%   the time/frequcency localization of the the analysis and the synthesis
%   window. The 'itersine' as a tight smooth window is can be considered as
%   a good tradeoff between the analysis and the sythesis window. Indeed,
%   both are equaly good. In this experiment, we are going to construct a
%   pair of dual window that are both more concentrated than the
%   'itersine'. For this experiment, the shift in time $a$ is set to $30$
%   and the number of frequency channel $M$ is set to $60$. To be compactly
%   supported, the itersine window has to be of length: $60$.
%
%   To find the pair of dual window, we first choose a well localized
%   window $g$ and we construct a dual by solving the following problem: 
%
%   .. gd  = argmin_x   gamma || nabla F x||_2^2 + mu || nabla x ||_2^2
%
%   ..     such that    x is a dual windows of g and is compaclty supported
%
%   .. math:: \begin{split}  \text{gd}  =  \text{arg} \min_x   & \gamma \|Fx\|_2^2 + \mu \| \nabla x \|_2^2  \\     \text{such that } & x \text{ is a dual windows of }g  \\     \text{and } & x \text{ is compaclty supported } \end{split}
%
%   with 
%         - $\mu$   :  Smoothing parameter
%         - $\gamma$:  Smoothing parameter in frequency
%         - $x$     :  Optimization variable
%         - $F$     :  The Fourier matrix 
%
%   Experiments
%   -----------
%
%   In this the experiment, we construct a pair of compactly supported dual
%   windows $g,gd$ such that $G(g,30,60)$ and $G(gd,30,60)$ form dual
%   frames. To improve over the itersine, we allow some freedom for the
%   support. We choose for $g$ a Nuttall window of length $L_g = 120$. This
%   window is slightly broader than the Itersine window $gi$ in 
%   time, but very good decay and localization in frequency, see the
%   figures above. For the dual window, we consider a support of length $360$.     
%
%   The dual window we obtain (using the regularization parameters 
%   $\mu_1 = 0.1$, $\mu_2 = 1$), shows roughly the same time localization 
%   as the Nuttall window or the tight window obtained before and is quite
%   similar to the latter in frequency. 
%
%   Results
%   -------      
%
%   .. figure::
%
%      Windows in time domain
%
%    
%
%   .. figure::
%
%      Windows in frequency domain
%
%
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
    paramplot.pathfigure = 'figures/beat_itersine/'; 
else
    paramplot.pathfigure = 'figures_color/beat_itersine/'; 
end
paramplot.position = [100 100 600 400];
paramplot.titleweight ='bold';
paramplot.save = GLOBAL_save; % save figures
paramplot.baw = baw; % plot everything in black and white



%%
% parameters of frame and length
a=30;
M=60;
Ldual=360;

% For testing the reconstruction
Llong=dgtlength(10*Ldual,a,M);
Lplot=Ldual;

window1='nuttall';
% windows initial point
g=firwin(window1,2*M);

gi=firwin('itersine',M);

%parameter for convex optimization
% those alpha should be smaller than the value of the windows !
alpha1=0.00; % weight of the L1 norm in time
alpha2=0.00; % weight of the L1 norm in frequency
mu=1; % smoothing parameter
gamma=0.1;
omega=0;
maxit=500; % maximum number of iteration
tol=1e-7; % tolerance to stop iterating

%%
fprintf('Solve the optimization problem \n');

% call optimization routine
gd=gabfirdual(Ldual,g,a,M,'alpha',alpha1,'beta',alpha2,'mu',mu,...
    'gamma',gamma,'omega',omega,'maxit',maxit,'tol',tol,'quiet');



gd=real(gd);


%%

if figpaper
    paramplot.position = [100 100 300 200];
    
    figure();
    plot_time(gi);
    title('(a)');
    save_name='iter_time';
    plotfig(save_name,paramplot );
    figure();
    magresp(gi,'fir','dynrange',dr,'1');
    title('(b)');
    setplotfreq();
    save_name='iter_freq';
    plotfig(save_name,paramplot );
   
    figure()
    plot_ambiguity(gi,dr);
    title('(c)')
    save_name='iter_amb';
    plotfig(save_name,paramplot );
    
    figure();
    plot_time(g);
    title('(d)');
    save_name='ana_time';
    plotfig(save_name,paramplot );
    figure();
    magresp(g,'fir','dynrange',dr,'1');
    title('(e)');
    setplotfreq();
    save_name='ana_freq';
    plotfig(save_name,paramplot );
   
    figure()
    plot_ambiguity(g,dr);
    title('(f)')
    save_name='ana_amb';
    plotfig(save_name,paramplot );
    
    
    figure();
    plot_time(gd);
    title('(g)');
    save_name='syn_time';
    plotfig(save_name,paramplot );
    figure();
    magresp(gd,'fir','dynrange',dr,'1');
    title('(h)');
    setplotfreq();
    save_name='syn_freq';
    plotfig(save_name,paramplot );  
    figure()
    plot_ambiguity(gd,dr);
    title('(i)')
    save_name='syn_amb';
    plotfig(save_name,paramplot );
    
else

    figure();

    hold off
    plot_time(fir2long(g,Lplot),'r');
    hold on;
    plot_time(fir2long(gi,Lplot),'b');
    plot_time(fir2long(gd,Lplot)/max(abs(gd)),'g');
    if baw
        set_baw_color();
    end
    legend(window1,'Itersine','Dual window','location', 'NorthEast')
    
    title('(a)')
    save_name='time';
    paramplot.baw=0;
    plotfig(save_name,paramplot );
    paramplot.baw=baw;
    figure()
    hold off
    magresp(fir2long(g,Lplot),'fir','dynrange',dr,'1','opts',{'r'});
    hold on;
    magresp(fir2long(gi,Lplot),'fir','dynrange',dr,'1','opts',{'b'});

    magresp(gd,'fir','dynrange',dr,'1','opts',{'g'});
    if baw
        set_baw_color();
    end
    title('(b)')
    legend(window1,'Itersine','Dual window')
    setplotfreq();
    save_name='freq';
    paramplot.baw=0;
    plotfig(save_name,paramplot );
    paramplot.baw=baw;
end

% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',gabdualnorm(g,gd,a,M,Llong));


% figure()
% plot(20*log10(abs(fftshift(fir2long(gd,Lplot)))+eps));
% title('Solution window - time domain - Log')


%%
Lcomp=dgtlength(2*M+Ldual+1,a,M);

[A,B] = gabframebounds(g,a,M)

fac = 4;

[crit_mat(1,:),crit_mat2(1,:)] = compute_criteria(gi/sqrt(A),Lcomp,fac,Ldual);
[crit_mat(2,:),crit_mat2(2,:)] = compute_criteria(g/sqrt(A),Lcomp,fac,Ldual);
[crit_mat(3,:),crit_mat2(3,:)] = compute_criteria(gd*sqrt(A),Lcomp,fac,Ldual);
crit_mat = crit_mat(:,[2,3,8,9,10,11,4,5,1,6,7]);

%%

XX = [100,100,1,1,10,10,1,1,10,1,1]'*[1,1,1];
crit = crit_mat.*XX';

XX2 = [100,100,1,1,1]'*[1,1,1];
crit2 = [crit_mat2.*XX2',crit(:,[1,2])];

%mat2tex(crit,'table')
mat2tex(crit2,'table')
