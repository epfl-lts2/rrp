function [sol, infos] = gsp_tv_deconvolution_noise(G, x0, h, sigma, param)

if nargin<5
    param = struct;
end

if ~isfield(param,'verbose'), param.verbose = 1; end


% Fidelity term for Tikonov an TV
paramproj.epsilon = sqrt(G.N)*sigma;
paramproj.verbose = param.verbose - 1;
paramproj.tight = 0;
ffid.prox = @(x,T) gsp_proj_b2_filterbank(x, T, G, h, x0, paramproj);
ffid.eval = @(x) eps;

% TV regularizer
paramtv.verbose = param.verbose -1;
ftv.prox = @(x,T) gsp_prox_tv(x,T,G,paramtv);
ftv.eval = @(x) gsp_norm_tv(G,x);

% Solve the problem
[sol, infos] = solvep(x0,{ffid,ftv}, param);


end