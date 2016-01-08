function [sol, infos] = gsp_wiener_inpainting(G,x0, A, At, psd, psd_noise, param)

if nargin<7
    param = struct;
end

if ~isfield(param,'nu'), param.nu = 1; end


% Fidelity term for Wiener optimization
ffid.grad = @(x) 2*At(A(x)-x0);
ffid.eval = @(x) norm(A(x)-x0,'fro')^2;
ffid.beta = 2*param.nu;

[sol, infos] = gsp_wiener_optimization(G, x0, ffid, psd, psd_noise, param);

end