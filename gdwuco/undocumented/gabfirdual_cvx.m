function gd=gabfirdual_cvx(Ldual,g,a,M,alpha1,alpha2,mu,gamma)


%length
Lg=length(g);
L=dgtlength(Lg+Ldual+1,a,M);
b=L/M;

%initial point
    gsmall=long2fir(g,M);
    gdsmall=gabdual(gsmall,a,M);

%Projection
    glong=fir2long(g,L);
    G=tfmat('dgt',glong,M,a);
    d=[a/M;zeros(a*b-1,1)];
    
% Constraint
    Lz=Ldual/2;
    ind=Lz:L-Lz;

% Matrix
    F=tfmat('fourier',L) ;
    Grad=eye(L)-circshift(eye(L),[0,-1]);
    
    

cvx_begin
    variable X(L)
    cvx_precision best
    minimize( alpha1*norm(X,1) +alpha2 *norm(F*X,1) + mu* norm(Grad*X,2) + gamma * norm(Grad*F*X,2))
    subject to
        G'*X == d
        X(ind) == 0
cvx_end


gd=long2fir(X,Ldual);
end

    
