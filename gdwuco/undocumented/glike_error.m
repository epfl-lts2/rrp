function [ d ] = glike_error( x,glike )
%GLIKE_ERROR Compute the distance the set glike

px = glike_proj(x,glike);

d = norm(x-px);


end

