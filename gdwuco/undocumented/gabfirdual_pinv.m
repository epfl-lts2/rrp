function [gd,Gcut]=gabfirdual_pinv(Ldual,g,a,M)

% Determine the window. The window /must/ be an FIR window, so it is
% perfectly legal to specify L=[] when calling gabwin
[g,info]=gabwin(g,a,M,[]);

% Determine L. L must be longer than L+Ldual+1 to make sure that no convolutions are periodic
L=dgtlength(info.gl+Ldual+1,a,M);
b=L/M;



% Determine an L that is so large that the system is underdetermined,
% even when reducing the number of variables
% Ldual must be less than a*b
% => b must be larger than ceil(Ldual/a)

Lfirst=ceil(Ldual/2);
Llast=Ldual-Lfirst;
glong=fir2long(g,L);

F=frame('dgt',glong,M,a);
G=frsynmatrix(F,L);




Gcut=G([1:Lfirst,L-Llast+1:L],:);
%size(Gcut)
gd=pinv(Gcut.')*[a/M;zeros(a*b-1,1)];
