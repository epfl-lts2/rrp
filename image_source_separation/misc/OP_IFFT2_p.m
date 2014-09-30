function [ F ] = OP_IFFT2_p( f )
%OP_IFFT2_p Compute the 2D ifft of f row by row.
%   First: resize squarely the row, then compute the 2D IFFT

[I,J]=size(f);

s_i=sqrt(I);
F=reshape(ifft2(reshape(f,s_i,s_i,J)),I,J);


end