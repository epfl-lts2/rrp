function bound = sp_spectral_frame_bound(G,g,p)


F = gsp_filterbank_matrix(G,g);
[A,B] = gsp_filterbank_bounds(G,g);
g12max = max(sqrt(sum(F.^2,1)));

bound = B^(min(1/2,1/p))/A^(max(1/2,1/p)) * g12max^(abs(1-2/p));

end