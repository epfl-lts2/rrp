%% In this file, we test the double projection
clear all
close all



amax=10;


% windows
window_type = 'rand'; % 'rand' or any ltfat windows
window_type = 'nuttall';


Matp=zeros(amax);
Matpf=zeros(amax);
Matpa=zeros(amax);
A=zeros(amax);
Mf=zeros(amax);
ML=zeros(amax);
MLg=zeros(amax);

for a=1:amax
    for M=a:amax
        Lg=dgtlength(1,a,M);
        if strcmp(window_type,'rand')
            g=rand(Lg,1); 
        else
            g=firwin(window_type,Lg);
        end
        Fal1=frame('dgt',g,a,M);
        [Bm,~]=framebounds(Fal1);
        if Bm < 1e-6
            fprintf('a is: %i    M is: %i    L is: %i   : not a Frame with the current window\n',a,M,L);
        else
            fprintf('a is: %i    M is: %i    Lg is: %i   ',a,M,Lg);
            for Ldual=1:a*Lg
                L=dgtlength(Lg+Ldual+1,a,M);
                b=L/M;

                glong=fir2long(g,L);
                Fal=frame('dgt',glong,M,a);
                G=frsynmatrix(Fal,L);
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
                    prox= @(x,T) fir2long(proj_dual(long2fir(x,Ldual),T,param_proj2),L); % set the prox

                    p=prox(x,0);
                catch
                   p=zeros(size(x)); 
                end


                error_p=norm(G'*p-d);
                if error_p<10e-10
                    Matp(a,M)=Ldual/L;
                    Matpf(a,M)=Ldual;
                    Matpa(a,M)=Ldual/a;
                    A(a,M)=a;
                    Mf(a,M)=M;
                    ML(a,M)=L;
                    MLg(a,M)=Lg;
                    fprintf('Ldual is: %i\n',Ldual);
                    break
                end
                if Ldual==a*Lg
                    fprintf('Ldual not found: %g\n',error_p);
                end
            end
        end
    end    
end
%%
imagesc(Matp)


