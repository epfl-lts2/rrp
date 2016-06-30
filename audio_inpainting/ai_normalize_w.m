function [win_n] = ai_normalize_w(win,hop)
%AI_NORMALIZE_W Normalizes the window so that it yields a tight STFT
% 
%   Input parameters:
%       win         : input window
%       hop         : hop-size
%   Ouput parameters:
%       win_n       : normalized window
%
%

% Authors: Nathanael Perraudin, Nicki Hollighaus
% Date   : June 2016

N = length(win);
K = floor(N/hop);
win_n = win.^2;
z = win_n;
for n = 1:K,%shift to the left
    z(1:end-n*hop) = z(1:end-n*hop) + win_n(n*hop+1:end);
end
for n = 1:K,%shift to the right
    z(n*hop+1:end) = z(n*hop+1:end) + win_n(1:end-n*hop);
end
win_n = win./sqrt(z);

end