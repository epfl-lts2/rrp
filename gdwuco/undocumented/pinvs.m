function X = pinvs(A,varargin)
%PINVS   Pseudoinverse for symetric matrices.
%   X = PINVS(A) compute the pseudo inverse of a symetric matrix
%


if isempty(A)     % quick return
  X = zeros(size(A'),class(A));  
  return  
end

% Be sure that the matrix is symetric
A = (A +A')/2;

n = size(A,1);


[Q,S] = eig(A);

Q = fliplr(Q);
   
s = flipud(real(diag(S)));

if nargin == 2
  tol = varargin{1};
else
  tol = n * eps(max(s));
end
r = sum(s > tol);
if (r == 0)
  X = zeros(size(A),class(A));
else
  s = diag(ones(r,1)./s(1:r));
  X = Q(:,1:r)*s*Q(:,1:r)';
end

end
