function [  ] = plotfftreal( ffts,fs,st )
%PLOTFFT Plot the fft in a nice way
%   Usage: plotfftreal( ffts, fs );
%          plotfftreal( ffts, fs, st );
%
%   Input arguments:
%       ffts    : fft of the signal s
%       fs      : sampling frequency
%       st      : stem (set it to 1 form stem plot -  default 0)
%   Ouput arguments:
%       none
%  

% Author: Nathanael Perraudin
% Date: 12 Mai 2014

if nargin<3
    st = 0;
end

N = length(ffts);

w = linspace(0,fs-fs/N,N);

if st
    stem(w,abs(ffts(:)));
else
    plot(w,real(ffts(:)));
end
xlabel('Frequency')
ylabel('Modulus')
drawnow;


end

