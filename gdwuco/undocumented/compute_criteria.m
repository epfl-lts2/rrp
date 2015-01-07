function [l2,gradl2,Fgradl2,s0,s0w,l1,Fl1,var,Fvar,Evar,FEvar, ...
    TW,ML,TWMLrat,SLatt,SLdec] = compute_criteria(g,Lcomp,fac,Ldual)

if nargin < 4
    Ldual = length(g);
if nargin < 3
    fac = 1;
    if nargin < 2
        Lcomp = length(g);
    end
end
end
%Lg = length(g);
LLg = fac*Lcomp;
g = fir2long(g,LLg);
Fg = fft(g)/sqrt(LLg);

l2 = norm(g,2);
l1 = norm(g,1);
Fl1 = norm(Fg,1)/sqrt(fac);

gradl2 = norm(g-circshift(g,1),2);
Fgradl2 = norm(Fg-circshift(Fg,1),2)*fac;

if mod(LLg,2)
    w=[0:1:(LLg-1)/2,(LLg-1)/2:-1:1]';
else
    w=[0:1:LLg/2-1,LLg/2:-1:1]';
end
w=w/sqrt(Lcomp);

var = norm(g.*w.^2,1);
Fvar = norm(Fg.*w.^2,1)/fac^(5/2);

Evar = norm(g.*w,2);
FEvar = norm(Fg.*w,2)/fac;

gau = pgauss(LLg,1/fac)/sqrt(Lcomp);
G = dgt(g,gau,1,LLg);
s0 = norm(G(:),1)/fac;


if mod(LLg,2)
    w=[0:1:(LLg-1)/2,(LLg-1)/2:-1:1]';
else         
    w=[0:1:LLg/2-1,LLg/2:-1:1]';
end
w=w/sqrt(LLg);
W = repmat(w,1,LLg).^2/fac + fac*repmat(w',LLg,1).^2;   
W=log(1+W);

s0w = norm(W(:).*G(:),1)/fac;


SLFg = abs(Fg./max(Fg));
SLFgrad = SLFg-circshift(SLFg,1);

G = g./max(g);
TW = find(G(1:ceil(end/2)-1)>=.5,1,'last');
TW = (2*TW-1)/Ldual;

ML = find(SLFgrad > 0,2);

SLatt = -max(20*log10(SLFg(ML(2):end-ML(2)+1)));

ML = (2*ML(find(ML>1,1,'first'))-3)/LLg;

TWMLrat = (TW*ML).^-1;

try
ind = intersect(find(SLFgrad>0),find(SLFgrad<0)-1);
SLdec = min(-SLatt-20*log10(SLFg(ind(floor(end/2)-2:floor(end/2)))));
catch
    SLdec = nan;
end


if nargout < 3
    l2 = [l2,gradl2,Fgradl2,s0,s0w,l1,Fl1,var,Fvar,Evar,FEvar];
    gradl2 = [TW,ML,TWMLrat,SLatt,SLdec];
end
