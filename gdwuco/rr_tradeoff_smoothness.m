%RR_TRADEOFF_SMOOTHNESS Smoothness time frequency tradeoff
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
%   It is well known that no function can be arbitrarily concentrated in
%   both the time and the frequency domain simultaneously. When choosing a
%   dual window to a given Gabor frame the concentration is further limited
%   by the duality conditions, the shape and the quality of the given
%   frame. Oversimplified, a badly conditioned Gabor frame (with large
%   frame bound ratio $B/A$), admits only badly concentrated duals.
%   However, even if the canonical dual window is well concentrated
%   overall, applications might benefit from the improvement of time
%   concentration versus frequency concentration and vice-versa. To see
%   this, just recall that the TF shape of the synthesis window limits the
%   precision of TF processing.
%
%   The following experiment demonstrates the surprising flexibility that
%   the set of dual windows allows when choosing the appropriate TF
%   concentration trade-off. The system parameters are the same as in
%   Experiment 1: $L_g = 240$, $a = 50$, $M=300$ and dual window support
%   less or equal to $L_h = 600$. However, to provide a more diverse set of
%   examples, we exchanged the Tukey window for an Itersine
%   window. We selected the time and frequency gradient priors to control
%   the TF spread and optimize
%   
%   .. argmin_x  lambda_1 || nabla x ||_2^2 + lambda_2 || nabla x ||_2^2
%   
%   ..          such that x is dual and compactly supported
%
%   .. math: \mathop{\operatorname{arg~min}}\limits _{x \in \mathcal{C}_{\text{dual}}\cap \mathcal{C}_{\text{supp}}} \lambda_1 \|\nabla \mathcal{F} x\|_2^2 + \lambda_2\|\nabla x\|_2^2
%   
%   for varying positive $\lambda_1,\lambda_2$, therefore balancing the
%   both concentration measures against one another. Recall that 
%   $|| \nabla Fx ||_2^2$ leads to concentration in time, while 
%   $||\nabla x ||_2^2$ promotes concentration in frequency.   
%
%   The results, presented in the figure below, demonstrate the large
%   amount of flexibility when choosing the TF concentration trade-off. It
%   also shows that extreme demands on either time or frequency
%   concentration come at the cost of other properties. In this particular
%   example, time concentration comes at the cost of worse sidelobe
%   attenuation, while strong demands on frequency concentration inhibit
%   quick frequency decay. Despite this, all four solution windows behave
%   as expected and show reasonable to very good overall TF concentration.
%
%   .. figure::
%
%      Analysis window in time domain
%
%      The chosen windows is a 'itersine'
%
%   .. figure::
%
%      Analysis window in frequency domain
%
%      
%
%   .. figure::
%
%      Results of optimization in the time domain
%
%      
%
%   .. figure::
%
%      Results of optimization in the frequency domain
%
% 
%
%   .. figure::
%
%      Results of optimization: ambiguity function
%
%      
%


%% initialization
clear;
close all;
% clc


global GLOBAL_save
global GLOBAL_baw
global GLOBAL_figpaper



% plotting
baw=GLOBAL_baw; % plot everything in black and white
dr = 160;
figpaper = GLOBAL_figpaper;


if GLOBAL_baw
    paramplot.pathfigure = 'figures/tradeoff_smoothness/'; 
else
    paramplot.pathfigure = 'figures_color/tradeoff_smoothness/'; 
end
paramplot.position = [100 100 600 400];
paramplot.titleweight ='bold';
paramplot.save=GLOBAL_save;
paramplot.baw=baw;



%% parameters of frame and length
a=50;
M=300;
L=300;


Ldual = 600;
Ltest = 8*M;

% windows
%g=fir2long(firwin('rect',L/4),L);
% g = fir2long(fftshift(tukeywin(240,.2)),L);
% g=fir2long(firwin('hann',240),L);
g=fir2long(firwin('itersine',240),L);

%g=pgauss(L,1);

% g=randn(size(g));
% g=g/norm(g);


%% Canonical dual
% The windows is tight. So the canonical dual is the same as the original
% window
fprintf('Compute the canonical dual \n');
if Ldual>M
gcan=fir2long(gabdual(g,a,M),Ldual);
else
gcan=long2fir(gabdual(g,a,M,L),Ldual);
end



% Verification: this number should be close to 0.
fprintf('  Reconstruction error of the canonical dual %g \n',...
                                gabdualnorm(g,gcan,a,M,Ltest));


%%

maxit=75;   % maximum number of iteration
tol=1e-3;    % tolerance to stop iterating




%% Solves the problem for 1/10

%parameter for convex optimization
mu=10;          % smoothing parameter in time
gamma=100;      % smoothing parameter in frequency


fprintf('Solve the optimization problem : 1/10\n');
% call optimization routine
gd=gabfirdual(Ldual,g,a,M,'maxit',4*maxit,'tol',0.01*tol,...
'mu',mu,'gamma',gamma,'evolution',1);
gd_grad110=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,Ltest));


%% Solves the problem for 5/1

%parameter for convex optimization
mu=500;         % smoothing parameter in time
gamma=100;      % smoothing parameter in frequency


fprintf('Solve the optimization problem : 3/1\n');
% call optimization routine
gd=gabfirdual(Ldual,g,a,M,'maxit',4*maxit,'tol',0.01*tol,...
'mu',mu,'gamma',gamma,'evolution',1);
gd_grad31=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,Ltest));
                            
%% Solves the problem for 20/1

%parameter for convex optimization
mu=1000;         % smoothing parameter in time
gamma=10;      % smoothing parameter in frequency


fprintf('Solve the optimization problem : 100/1\n');
% call optimization routine
gd=gabfirdual(Ldual,g,a,M,'maxit',4*maxit,'tol',0.01*tol,...
'mu',mu,'gamma',gamma,'evolution',1);
gd_grad201=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,Ltest));
                            
%% Solves the problem for 1000/1

%parameter for convex optimization
mu=1000;         % smoothing parameter in time
gamma=1;      % smoothing parameter in frequency


fprintf('Solve the optimization problem : 1000/1\n');
% call optimization routine
gd=gabfirdual(Ldual,g,a,M,'maxit',4*maxit,'tol',0.01*tol,...
'mu',mu,'gamma',gamma,'evolution',1);
gd_grad10001=real(gd);
% Verification: this number should be close to 0.
fprintf('  Reconstruction error of optimization problem: %g \n',...
                                gabdualnorm(g,real(gd),a,M,Ltest));

                            



%% Plotting the results

if figpaper
    paramplot.position = [100 100 200 165];

    figure();
    plot_time(g);
    title('(a)');
    save_name='analysis_time';
    plotfig(save_name,paramplot );
    figure()
    magresp(g,'fir','dynrange',dr,'1');
    title('(b)');
    setplotfreq();
    save_name='analysis_freq';
    plotfig(save_name,paramplot );


    figure()
    plot_time(gd_grad10001);
    title('(a)');
    ylatex('1000 ||\nabla x||_2^2 + ||\nabla F x||_2^2');
    save_name='gd_10001_time';
    plotfig(save_name,paramplot );
    
    figure()
    plot_time(gd_grad201);
    title('(d)');
    ylatex(' 100 ||\nabla x||_2^2 + ||\nabla F x||_2^2');
    save_name='gd_201_time';
    plotfig(save_name,paramplot );
   
    figure
    plot_time(gd_grad31);
    title('(g)');
    ylatex(' 5 ||\nabla x||_2^2 + ||\nabla F x||_2^2');
    save_name='gd_31_time';
    plotfig(save_name,paramplot );
    
    figure
    plot_time(gd_grad110);
    title('(j)');
    ylatex(' ||\nabla x||_2^2 + 10 ||\nabla F x||_2^2');
    save_name='gd_110_time';
    plotfig(save_name,paramplot );




    figure()
    magresp(gd_grad10001,'fir','dynrange',dr,'1');
    title('(b)');
    setplotfreq();
    save_name='gd_10001_freq';
    plotfig(save_name,paramplot );
    figure()
    magresp(gd_grad201,'fir','dynrange',dr,'1');
    title('(e)');
    setplotfreq();
    save_name='gd_201_freq';
    plotfig(save_name,paramplot );
    figure()
    magresp(gd_grad31,'fir','dynrange',dr,'1');
    title('(h)');
    setplotfreq();
    save_name='gd_31_freq';
    plotfig(save_name,paramplot );
    figure()
    magresp(gd_grad110,'fir','dynrange',dr,'1');
    title('(k)');
    setplotfreq();
    save_name='gd_110_freq';
    plotfig(save_name,paramplot );


    figure()
    plot_ambiguity(gd_grad10001,dr);
    title('(c)')
    save_name='ambiguity_10001';
    plotfig(save_name,paramplot );
    
    figure()
    plot_ambiguity(gd_grad201,dr);
    title('(f)');
    save_name='ambiguity_201';
    plotfig(save_name,paramplot );  
    
    figure()
    plot_ambiguity(gd_grad31,dr);
    title('(i)');
    save_name='ambiguity_31';
    plotfig(save_name,paramplot );  
    
    figure()
    plot_ambiguity(gd_grad110,dr);
    title('(l)');
    save_name='ambiguity_110';
    plotfig(save_name,paramplot );

    
else
    paramplot.position = [100 100 300 200];
    figure(1);
    plot_time(g);
    title('(a)');
    save_name='analysis_time';
    plotfig(save_name,paramplot );
    figure(2)
    magresp(g,'fir','dynrange',dr,'1');
    title('(b)');
    setplotfreq();
    save_name='analysis_freq';
    plotfig(save_name,paramplot );

    figure(3);
    hold on;

    plot_time(gd_grad10001,'g');
    plot_time(gd_grad110,'k');
    plot_time(gd_grad31,'b');
    plot_time(gd_grad201,'r');


    if baw
        set_baw_color();
    end

    legend('(1000/1)','(1/10)','(5/1)',...
        '(100/1)');
    title('(c)');
    setplottime(gd_grad110);
    %xlim([0,300]);
    save_name='comp_opt_time';
    paramplot.baw=0;
    plotfig(save_name,paramplot );


    paramplot.baw=0;
    figure(4)
    hold on;
    magresp(gd_grad10001,'fir','dynrange',dr,'1','opts',{'g'});
    magresp(gd_grad110,'fir','dynrange',dr,'1','opts',{'k'});
    magresp(gd_grad31,'fir','dynrange',dr,'1','opts',{'b'});
    magresp(gd_grad201,'fir','dynrange',dr,'1','opts',{'r'});



    if baw
        set_baw_color();
    end


    legend('(1000/1)','(1/10)','(5/1)',...
        '(100/1)');
    title('(d)');
    setplotfreq();
    save_name='comp_opt_freq';
    plotfig(save_name,paramplot );

    paramplot.baw=baw;
    figure(5);
    paramplot.position = [100 100 600 600];
    save_name='ambiguity';

    subplot(221)
    plot_ambiguity(gd_grad110,dr);
    title('(e)');
    plotfig(save_name,paramplot );

    subplot(222)
    plot_ambiguity(gd_grad31,dr);
    title('(f)');
    plotfig(save_name,paramplot );

    subplot(223)
    plot_ambiguity(gd_grad201,dr);
    title('(g)');
    plotfig(save_name,paramplot );

    subplot(224)
    plot_ambiguity(gd_grad10001,dr);
    title('(h)')
    plotfig(save_name,paramplot );
end
    

%%
Lcomp=dgtlength(L+Ldual+1,a,M);

fac = 4;
                            
[crit_mat(1,:),crit_mat2(1,:)] = compute_criteria(gcan,Lcomp,fac);
[crit_mat(2,:),crit_mat2(2,:)] = compute_criteria(gd_grad110,Lcomp,fac);
[crit_mat(3,:),crit_mat2(3,:)] = compute_criteria(gd_grad31,Lcomp,fac);
[crit_mat(4,:),crit_mat2(4,:)] = compute_criteria(gd_grad201,Lcomp,fac);
[crit_mat(5,:),crit_mat2(5,:)] = compute_criteria(gd_grad10001,Lcomp,fac);
% Y = (crit_mat - repmat(min(crit_mat(1:5,:)),5,1)); 
% Y = Y./repmat(max(Y(1:5,:)),5,1);
crit_mat = crit_mat(:,[2,3,8,9,10,11,4,5,1,6,7]);

%Y
%%

XX = [1000,1000,10,10,100,100,1,1,100,10,10]'*[1,1,1,1,1];
crit = crit_mat.*XX'

XX2 = [100,100,1,1,1]'*[1,1,1,1,1];
crit2 = [crit_mat2.*XX2',crit(:,[1,2])]     

mat2tex(crit([2,3,4,5],:),'table')
mat2tex(crit2([2,3,4,5],:),'table')

                            
