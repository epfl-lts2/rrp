function [ g ] = aritight( a,M )
%ARITIGHT Create the aritight windows for parameters a and M


Lg=M;
Ldual=M;

% For testing the reconstruction
Llong=dgtlength(10*(Ldual+Lg),a,M);


% Initial windows
g=firwin('itersine',Lg);


%parameter for convex optimization
    mu=0;         % smoothing parameter in time
    gamma=0;      % smoothing parameter in frequency
    delta=0.1/sqrt(2*M);  % parameter for the Gabor transform



maxit=500;   % maximum number of iteration
tol=1e-7;    % tolerance to stop iterating


fprintf('Compute the window \n');
% call optimization routine
gd=gabfirdual(Ldual,g,a,M,'mu',mu,...
    'gamma',gamma,'maxit',maxit,'tol',tol,'delta',delta);




%%
g=real(fir2long(gabtight(gd,a,M),Ldual));

if gabdualnorm(g,g,a,M,Llong)>10e-14
fprintf('  Reconstruction error of the tight window: %g \n',...
                                gabdualnorm(g,g,a,M,Llong));
end


end

