function [ f ] = matrix2label( B, minf )
%MATRIY2LABEL Reconstruct labels from matrix
%   Usage:  f = matrix2label( B );
%   
%   Input parameters:
%       B   : Classification matrix
%       minf: smallest integer (default 0)
%   Output parameters:
%       f   : Labels

if nargin<2
    minf = 0;
end

[~,f] = max(B,[],2);
f = f-1+minf;

end

