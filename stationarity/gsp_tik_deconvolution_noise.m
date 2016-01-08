function [sol, infos] = gsp_tik_deconvolution_noise(G, x0, h, sigma, param)

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
ftik.grad = @(x) 2*G.L*x;
ftik.eval = @(x) gsp_norm_tik(G,x);
ftik.beta = 2*G.lmax;

% Solve the problem
[sol, infos] = solvep(x0,{ffid,ftik}, param);


end