function p = simplex_proj(s,param)
%projection on the simplex being the intersection of the unit diamond and
%the positiv orthant.
%s : of size NxI with N: nb samples, I: nb dimensions
%p : the projection of s.
%param is optional
%param.algo = {'duchi','dykstra'}
%param.err for Dykstra only
%
%simon Arberet 08/02/2011

[N,I] = size(s);

if nargin<2
    param.algo = 'duchi';
    %algorithm of section 3 of John Duchi et al.
    %"Efficient Projections onto the ?1-Ball for Learning in High
    %Dimensions",ICML,2008.
end

if strcmp(param.algo,'dykstra')
    if isfield(param,'err')
      err = param.err;
    else
      err = 1e-6;
    end
    p = s;
    p_prev = p;
    %s1 = s;
    z = zeros(N,I);
    z1 =z;
    iter11 = 1;
    err11 = 1;
    
    %projection to convex constraint
    while (err11>err) && (iter11<100)
        
        %Dykstra's algorithm to project onto the intersection of convex
        %sets
        t = p - z1;
        s1 = t;
        s1(s1<0) = 0;
        z1 = s1 - t;
        
        t = s1 - z;
        SUM = sum(t,2);
        d = (1-SUM)/I; % L1 = 1
        p = t + repmat(d,1,I);
        z = p - t;
        
        err11 = norm(p_prev - p, 'fro')/norm(p_prev,'fro');
        p_prev = p;
        iter11 = iter11+1;
    end
    sprintf('Dykstra algorithm: %d iterations',iter11-1)
else %exact
    z=1;
    Mu = sort(s',1,'descend');
    Thetas = repmat((1./(1:I))',[1 N]).*(cumsum(Mu)-z);
    W = Mu - Thetas > 0;
    rho = sum(W);
    theta = Thetas(sub2ind([I,N],rho,1:N));
    w = max(s' - repmat(theta,[I 1]),0);
    p=w';
end