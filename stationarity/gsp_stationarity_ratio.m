function r = gsp_stationarity_ratio(G, C)
%GSP_STATIONARITY_RATIO Assert the stationarity level of some data
%   Usage:  r = gsp_stationarity_ratio(G, C)
%
%   Input parameters:
%         G          : Graph
%         C          : Covariance matrix
%   Output parameters:
%         r          : Ratio
%
%   This function compute the ratio of energy contained into the diagonal
%   of the Fourier covariance matrix:
%
%   .. T = U' * C * U 
%
%   .. math:: \Tau = U^{*} C U
%
%   References: perraudin2016stationary


% Author : Nathanael Perraudin
% Date: 6 January 2016

CF = G.U' * C * G.U;

r = gsp_diagonal_ratio(CF);



end


function r = gsp_diagonal_ratio(M)

if size(M,1) == 1 || size(M,2) ==1
    error('This function acts on martrices.')
end

r = norm(diag(M))/norm(M,'fro');

end