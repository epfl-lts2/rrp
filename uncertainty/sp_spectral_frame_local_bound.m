function [bound, approx_bound, ktilde, itilde, dd] = sp_spectral_frame_local_bound(G,g,p,ii,k)

    bound = zeros(size(k));
    approx_bound = zeros(size(k));
    ktilde = zeros(size(k));
    itilde = zeros(size(k));
    dd = zeros(size(k));


    [A,B] = gsp_filterbank_bounds(G,g);
    Nf = numel(g);
    for jj = 1:numel(k)
        gi = gsp_localize(G,g{k(jj)},ii)/sqrt(G.N);

        Agg = gsp_filter_analysis(G,g,gi);

        [mA,k2] = max(abs(gsp_vec2mat(Agg,Nf)'));
        [~, itilde(jj)] = max(mA);
        ktilde(jj) = k2(itilde(jj));
        dd(jj) = gsp_hop_distanz(G,ii,itilde(jj));
        approx_bound(jj) = (B*G.N)^(1/p)/sqrt(A*G.N)*(norm(gi,2))^(abs(1-2/p));      
        gki = gsp_localize(G,g{ktilde(jj)},itilde(jj))/sqrt(G.N);
        bound(jj) = (B)^(1/p)/sqrt(A)*(norm(gki,2))^(abs(1-2/p)); 
    end
    
end