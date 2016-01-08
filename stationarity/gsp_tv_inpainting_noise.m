function [sol, infos] = gsp_tv_inpainting_noise(G, x0, Mask, sigma, param)

if nargin<5
    param = struct;
end

if ~isfield(param,'verbose'), param.verbose = 1; end


% Fidelity term for Tikonov an TV
paramproj.A = @(x) Mask.*x;
paramproj.At = @(x) Mask.*x;
paramproj.y = x0;
paramproj.epsilon = sqrt(sum(Mask(:)))*sigma;
paramproj.verbose = param.verbose - 1;
paramproj.tight = 1;
ffid_classic.prox = @(x,T) proj_b2(x,T,paramproj);
ffid_classic.eval = @(x) eps;

% TV regularizer
paramtv.verbose = param.verbose -1;
ftv.prox = @(x,T) gsp_prox_tv(x,T,G,paramtv);
ftv.eval = @(x) gsp_norm_tv(G,x);

% Solve the problem
[sol, infos] = solvep(x0,{ffid_classic,ftv}, param);


end