function Theta = theta (N)
% theta - CREATE A RANDOM +,-1 VECTOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theta = theta (N) creates a random +,-1 vector for random convolution
% Compressed Sensing
% with NB_MEAS 1.
% 
% Author: Mahdad Hosseini Kamal
% Contact: mahdad.hosseinikamal@epfl.ch
% Date: May 2010
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Theta = rand(N,1);
Theta(Theta>0.5) = 1;
Theta(Theta<0.5) = -1;
end
