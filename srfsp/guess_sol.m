function [ guess ] = guess_sol( sol )
%GUESS_SOL Guess the solution

Ntot = length(sol);

indbig = find(abs(sol) > 0.1*abs(mean(sol))*sqrt(Ntot));

guess = zeros(Ntot,1);


curr_ind = indbig(1);

for ii=2:length(indbig)
    if (indbig(ii) ~= indbig(ii-1)+1) || (ii == length(indbig))
        localpic = sol(curr_ind:indbig(ii-1));
        ng = norm(localpic,1);
        [~,ind_local] = max(abs(localpic));
        guess(curr_ind+ind_local-1) = ng;        
        
        curr_ind = indbig(ii);
    end
    
    
end

guess = guess.* exp(1i*angle(sol));

end

