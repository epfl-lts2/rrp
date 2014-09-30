%% Minimal support experiment.
%
%   This file computes the minimum support for a given $a$ and $M$ parameter.
%   The length of the analysis window is computed directly from a and $M$ as
%   the smaller valid length for those parameters::
%
%               L=dgtlength(1,a,M);
%
%   For given window, we then compute the WR equation system and we force
%   short support to the dual windows by truncation. 
%   
%   * To be sure that a dual windows exist, we compute the projection of
%     one random signal and check that there is no reconstruction error
%
%   * To find the minimal support, we start with a length of 1 and increase
%     the support until we have perfect reconstruction.
%
%   * At the end of the experiment, I compute the number of equation in the
%     WR system and compare it with the minimal length of the dual window.
%
%   Personnally I do not exactly get what happen in this experiment. Most
%   of the Gabor transform that I construct are not frame but I can still
%   get a very good reconstruction. I guess we cannot chose $a$ and $M$
%   independently.
%
%   
%

%clear all
close all


amax=30;
amin=1;

window_type = 'rand'; % 'rand' or any ltfat windows
window_type = 'nuttall';


% windows


Matp=zeros(amax);
Matpf=zeros(amax);
Matpa=zeros(amax);
A=zeros(amax);
Mf=zeros(amax);
ML=zeros(amax);

for a=amin:amax
    
    for M=a:amax
        L=dgtlength(1,a,M);
        b=L/M;
        if strcmp(window_type,'rand')
            g=rand(L,1); 
        else
            g=firwin(window_type,L);
        end
        glong=fir2long(g,L);
        Fal=frame('dgt',glong,M,a);
        Fal1=frame('dgt',glong,a,M);
        [Bm,~]=framebounds(Fal1);
        if Bm < 1e-7
            fprintf('a is: %i    M is: %i    L is: %i   : not a Frame with the current window\n',a,M,L);
        else
            G=frsynmatrix(Fal,L);
            fprintf('a is: %i    M is: %i    L is: %i   ',a,M,L);
            for Ldual=1:L
                d=[a/M;zeros(a*b-1,1)];
                Lfirst=ceil(Ldual/2);
                Llast=Ldual-Lfirst;
                Gcut=G([1:Lfirst,L-Llast+1:L],:);

                x=glong/norm(glong);

                param_proj2.verbose=0;
                param_proj2.y=d;
                param_proj2.A=Gcut';
                try
                    param_proj2.AAtinv=pinvs(Gcut'*Gcut);
                    prox= @(x,T) fir2long(proj_dual(long2fir(x,Ldual),T,param_proj2),L);            
                    p=prox(x,0);
                catch
                    p=zeros(size(x));
                end
                error_p=norm(G'*p-d);
                if error_p<10e-8
                    Matp(a,M)=Ldual/L;
                    Matpf(a,M)=Ldual;
                    Matpa(a,M)=Ldual/a;
                    A(a,M)=a;
                    Mf(a,M)=M;
                    ML(a,M)=L;
                    fprintf('Ldual is: %i\n',Ldual);
                   % if floor(Ldual/a)==Ldual/a
                        break
                    %end
                end
                if Ldual==L
                    fprintf('Ldual not found: %g\n',error_p); 
                    gd = gabdual(g,a,M,L);
                    gabdualnorm(g,gd,a,M,L)
                    norm(G'*gd-d)
                    error('Strange....')
                end
            end
        end
    end    
end
%%
imagesc(Matp)
figure();

imagesc(Matpa)

ML./Mf.*A-Matpf %   ML./Mf.*A == number of linear equation
                %   and it is also equal to the minimum dual lenght for
                %   some windows for other not.
