function sol = glike_proj(x,glike)
%GLIKE_PROJ project x onto the glike set


glike = fir2long(glike,length(x));

% A) compute the closest point c*glike to x

c = sum(glike.^2) ./ sum(glike.*x);

% B) Return the projection

sol = c*glike;

end