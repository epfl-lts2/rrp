function sol = grm_estimator(C, M, s, sigma)

if nargin<4
    sigma = 0;
end

if size(M,1)<numel(M)
    error('M should be a column vector')
end

Cnn = C(logical(M),logical(M));
%Cuu = C(~logical(M),~logical(M));
Cun = C(~logical(M),logical(M));
sol = zeros(size(s));
sol(logical(M),:) = Cnn * (( Cnn + sigma*eye(sum(double(M))) ) \ s(logical(M),:));
sol(~logical(M),:) = (Cun * (( Cnn + sigma*eye(sum(double(M))) ) \ s(logical(M),:)));

end