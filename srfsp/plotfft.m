function [  ] = plotfft( ffts,fs )
%PLOTFFT Plot the fft in a nice way
%   Usage: plotfft( ffts );
%
%   Input arguments:
%       ffts    : fft of the signal s
%       fs      : sampling frequency
%   Ouput arguments:
%       none
%  

% Author: Nathanael Perraudin
% Date: 12 Mai 2014

N = length(ffts);

w = linspace(0,fs-fs/N,N);

plot(w,abs(ffts(:)));
xlabel('Frequency')
ylabel('Modulus')
xlim([0,fs/2]);
drawnow;

end

