function [ F ] = OP_FFT2_p( f )
%OP_FFT2_p Compute the 2D fft of f row by row.
%   First: resize squarely the row, then compute the 2D FFT

[I,J]=size(f);

s_i=sqrt(I);
F=reshape(fft2(reshape(f,s_i,s_i,J)),I,J);


end

