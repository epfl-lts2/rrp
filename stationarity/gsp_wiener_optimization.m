function [sol, infos] = gsp_wiener_optimization(G, x0, ffid, psd, psd_noise, param)

if nargin<6
    param = struct;
end

if isnumeric(psd_noise)
    wl = @(x) psd_noise./(psd(x)+eps);
else
    wl = @(x) psd_noise(x)./(psd(x)+eps);
end

f = @(x) 1./(wl(x)+1);

% Wiener term 
fwiener.prox = @(x,T) gsp_filter_analysis(G,f,x, param);
fwiener.eval = @(x) 0.5*norm(gsp_filter_analysis(G,wl,x,param),'fro')^2;

% Call the solver
[sol , infos ] = solvep(x0,{ffid,fwiener},param);

end