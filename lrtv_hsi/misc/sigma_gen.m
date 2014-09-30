function W_vector = sigma_gen(n)
%CREATE RANDOM PHASE SHIFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M = sigma_gen(n)
%
%
% Genarate the Sigma matrix for Random Convolution
% n is the number of pixels inside the image
% M is a sparse diagonal matrix which its terms are random phase
%
%
%         --                        --
%         |sigma_1 0 ...             |
%         |0      sigma_2 ...        |
% SIGMA = |.   .                     |
%         |.    .                    |
%         |.     .                   |
%         |0                sigma_n  |
%         --                        --
%
%
% w = 1        --> sigma_1     +,-1 with equal prob
% 2<=w<n/2+1   --> sigma_w     exp(jtheta_w) theta_w Uniform([0,2pi])
% w = n/2+1    --> sigma_n/2+1 +,-1 with equal prob
% n/2+2<=w<=n  --> sigma_w     conjugate of sigma_n-w+2 
%
% because of huge amount of data we just store the diagonal
%
%
%            --       --
%            |sigma_1  |
% W_vector = |sigma_2  |
%            | .       |
%            | .       |
%            | .       |
%            |sigma_n  |
%            --       --
%
%
% Author: Mahdad Hosseini Kamal
% Contact: mahdad.hosseinikamal@epfl.ch
% Date: May 2010
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W_vector = zeros(n,1);

if(rand>.5)
    W_vector(1,1) = 1;
else
    W_vector(1,1) = -1;
end

% for ii=n/2:-1:2
%     W_vector(ii,1) = exp(1j*2*pi*rand(1));
%     W_vector(n-ii+2,1) = conj(W_vector(ii));
% end

 W_vector(2:n/2,1)   = exp(1j*2*pi*rand(n/2-1,1));
 W_vector(n/2+2:n,1) = conj(W_vector(n/2:-1:2,1));

if(rand>.5)
    W_vector(n/2+1,1) = 1;
else
    W_vector(n/2+1,1) = -1;
end

%SIGMA = sparse(diag(W_vector));

end