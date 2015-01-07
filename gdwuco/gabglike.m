function [gd,relres,iter] = gabglike(g, glike, a,M,varargin)
%GABGLIKE Compute a dual window using convex optimization
%   Usage: gd=gabglike(g,glike,a,M);
%          gd=gabglike(g,glike,a,M, varagin);
%
%   Input parameters:
%     g      : Window function /initial point (tight case)
%     glike  : Design function
%     a      : Time shift
%     M      : Number of Channels
%
%   Output parameters:
%     gd     : Dual window
%     iter   : Number of iteration
%     relres : Reconstruction error
%
%   `gabconvexopt(g,a,M)` computes an window *gd* which is the optimal
%   solution of the convex optimization problem below
%
%   .. gd  = argmin_x    ...
%
%   ..     such that  x satifies the constraints
%
%   .. math:: \begin{split}  \text{gd}  = & \text{arg} \min_x   & ... \\ &    \text{such that }& x \text{ satisfies the constraints} \end{split}
%
%   Three constraints are possible:
%   
%   * *x* is dual with respect of g
%
%   * *x* is tight
%
%   * *x* is compactly supported on Ldual
%
%   WARNING: this function require the unlocbox! You can download it at
%   unlocbox.sourceforge.net
%
%   The function uses an iterative algorithm to compute the approximate.
%   The algorithm can be controlled by the following flags:
%
%     'glike',g_l  $g_l$ is a windows in time. The algorithm try to shape
%                  the dual window like $g_l$. Normalization of $g_l$ is done
%                  automatically. To use option omega should be different
%                  from 0. By default $g_d=0$.
%
%     'support' Ldual  Add a constraint on the support. The windows should
%                  be compactly supported on Ldual.
%
%     'tight'      Look for a tight windows
%
%     'dual'       Look for a dual windows (default)
%
%     'painless'   Construct a starting guess using a painless-case
%                  approximation. This is the default
%
%     'zero'       Choose a starting guess of zero.
%
%     'rand'       Choose a random starting phase.
%
%     'tol',t      Stop if relative residual error is less than the 
%                  specified tolerance.  
%
%     'maxit',n    Do at most n iterations. default 10
%
%     'print'      Display the progress.
%
%     'debug'      Display all the progresses.
%
%     'quiet'      Don't print anything, this is the default.
%
%     'fast'       Fast algorithm, this is the default.
%
%     'slow'       Safer algorithm, you can try this if the fast algorithm
%                  is not working. Before using this, try to iterate more.
%   
%     'evolution',p  Plot the window every p iteration (0 no plot), default 0
%
%     'gif',filename Produce a gif with the evolution of the window over
%                  iterations and save it to filename.gif. This function
%                  require 'evolution',p with p>0.
%
%     'hardconstraint' Force the projection at the end (default)
%
%     'softconstaint' Do not force the projection at the end
%
%   See also: gaboptdual, gabdual, gabtight, gabfirtight, gabopttight
   


% Author: Nathanael Perraudin
% Date  : 3 April 2014


if nargin<4
  error('%s: Too few input parameters.',upper(mfilename));
end;

if numel(g)==1
  error('g must be a vector (you probably forgot to supply the window function as input parameter.)');
end;

definput.keyvals.L=[];
definput.keyvals.lt=[0 1];
definput.keyvals.tol=1e-8;
definput.keyvals.maxit=1000;
definput.flags.print={'quiet','print','debug'};
definput.flags.algo={'fast','slow'};
definput.flags.constraint={'hardconstraint','softconstaint'};
definput.flags.startphase={'painless','zero','rand'};
definput.flags.type={'dual','tight'};
definput.keyvals.support=0;

definput.keyvals.evolution=-1;
definput.keyvals.gif=0;


[flags,kv]=ltfatarghelper({'L','tol','maxit'},definput,varargin);

% Determine the window. The window /must/ be an FIR window, so it is
% perfectly legal to specify L=[] when calling gabwin
[g,info]=gabwin(g,a,M,[],kv.lt,'callfun',upper(mfilename));

if kv.support
    Ldual=kv.support;
    % Determine L. L must be longer than L+Ldual+1 to make sure that no convolutions are periodic
    L=dgtlength(info.gl+Ldual+1,a,M);
else
    L=length(g);
    Ldual=L;
end

b=L/M;

% Determine the initial guess
if flags.do_zero
  gd_initial=zeros(Ldual,1);
end;

if flags.do_rand
  gd_initial=rand(Ldual,1);
end;

if flags.do_painless
  gsmall=long2fir(g,M);
  gdsmall=gabdual(gsmall,a,M);
  if length(gdsmall) <= Ldual
      gd_initial=fir2long(gdsmall,Ldual);
  else
      gd_initial=long2fir(gdsmall,Ldual);
  end
end;


% -------- do the convex optimization stuff

% Define the long original window
glong=fir2long(g,L);




%gabframebounds(g,a,M)


% Initial point
xin=gd_initial;
xin=fir2long(xin,L);


% -- * Setting the different prox for ppxa *--
% ppxa will minimize all different proxes

% value test for the selection constraint
F = {};

% - DUAL OR TIGHT?-    
if flags.do_tight
    % tight windows
    g2.prox= @(x,T) gabtight(x,a,M); % set the prox
    g2.eval= @(x) norm(x-gabdual(x,a,M,L)); % objectiv function
    
else
% - projection on a B2 ball -
    % Frame-type matrix of the adjoint lattice
    %G=tfmat('dgt',glong,M,a);
    Fal=frame('dgt',glong,M,a);
%    G=framematrix(Fal,L);
    G=frsynmatrix(Fal,L);
    d=[a/M;zeros(a*b-1,1)];
    

    % Using a direct projection (better solution)
    param_proj.verbose=flags.do_debug;
    param_proj.y=d;
    param_proj.A=G';
    param_proj.AAtinv=(G'*G)^(-1);
    g2.prox= @(x,T) proj_dual(x,T,param_proj); % set the prox
    g2.eval= @(x) norm(G'*x-d); % objectiv function
end
    


glike=glike/norm(glike);

g8.prox= @(x,T) glike_proj(x,glike);
g8.eval= @(x) glike_error(g2.prox(x,0),glike); % the objectiv function is the l2 norm
F{end+1} = g8;


% SUPPORT CONSTRAINT
if kv.support


% - function apply the two projections thanks to a poc algorithm.
    if flags.do_tight
    % - set null coefficient    
        g4.prox = @(x,T) forceeven(fir2long(long2fir(x,Ldual),L));
        g4.eval = @(x) 0;
        
        G={g2,g4};
        paramPOCS.tol=20*eps;
        paramPOCS.maxit=5000;
        paramPOCS.verbose=flags.do_print+flags.do_debug;
        paramPOCS.abs_tol=1;
        g5.prox = @(x,T) pocs(x,G,paramPOCS);
        % g5.prox = @(x,T) ppxa(x,G,paramPOCS);
        % g5.prox = @(x,T) douglas_rachford(x,g2,g4,paramPOCS);
        % g5.prox = @(x,T) pocs2(x,g2,g4,20*eps,2000, flags.do_print+flags.do_debug);
        g5.eval = @(x) 0;
        
        if flags.do_fast
            F{end+1} = g4;    
            F{end+1} = g2;   
        else
            F{end+1} = g5;               
        end
    else
        Fal=frame('dgt',glong,M,a);
        %G=framematrix(Fal,L);
        G=frsynmatrix(Fal,L);
        d=[a/M;zeros(a*b-1,1)];
        Lfirst=ceil(Ldual/2);
        Llast=Ldual-Lfirst;
        Gcut=G([1:Lfirst,L-Llast+1:L],:);
        param_proj2.verbose=flags.do_debug;
        param_proj2.y=d;
        param_proj2.A=Gcut';
        try
%            param_proj2.AAtinv=pinv(Gcut'*Gcut);
            param_proj2.AAtinv=pinv(Gcut)*pinv(Gcut');           
        catch
            error('GABOPTDUAL: Cannot pinv the matrix... We are sorry this is a Matlab bug');
        end
        g5.prox= @(x,T) fir2long(proj_dual(long2fir(x,Ldual),T,param_proj2),L); % set the prox
        g5.eval= @(x) norm(G'*x-d); % objectiv function
        F{end+1} = g5;    
    end
else
    F{end+1} = g2;
    g5 = g2;
end
   




  


  
    
    
% -- * PPXA function, the solver * --


% parameter for the solver
    param_solver.maxit=kv.maxit; % maximum number of iteration
    param_solver.tol=kv.tol;
    
    if flags.do_quiet
        param_solver.verbose=0;
    elseif flags.do_debug
        param_solver.verbose = 2;
    elseif param_solver.do_print
        param_solver.verbose = 1;
    end
    
    
    % solving the problem
    

    if kv.evolution > 0
        fig=figure(100);
        param_solver.do_sol=@(x) ...
            plot_window_proj(x,ceil(kv.evolution),fig,...
            @(x) g5.prox(x,0),kv.gif);
    end
    [gd,iter,~]=pocs(xin,F,param_solver);

    % Force the hard constraint
    if flags.do_hardconstraint
         gd=g5.prox(gd,0);
    end

    if kv.evolution > 0
        try 
            close(fig);
        catch
            warning('GABCONVEXOPT: You have messed with the figures!');            
        end
    end

    
    
    % compute the error
    if flags.do_tight
        relres=gabdualnorm(gd,gd,a,M,L);
    else
        relres=gabdualnorm(g,gd,a,M,L);
    end
    

   if kv.support
        % set the good size
        gd=long2fir(gd,Ldual);
   end

end



function x=forceeven(x)
% this function force the signal to be even
   x=  (x+involute(x))/2;
end
