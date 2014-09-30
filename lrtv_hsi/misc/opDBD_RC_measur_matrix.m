function op = opDBD_RC_measur_matrix(x_size, nb_channel, nb_meas)

%
% OPDBD_RC_MEASUR_MATRIX creates an operator of Distict Block Diagonal Random convolution measurement
% matrix
%
%   x_size:         input size
%   nb_channel:     number of channels
%   nb_meas:        number of measures
%
% Paper:
%           Compressive sensing 
% Author: Mohammad Golbabaee
% Contact: mohammad.golbabaei@epfl.ch
% Date: Oct. 2011


 W_vector = zeros(x_size, nb_channel);
 for i = 1:nb_channel
 W_vector(:,i) = sigma_gen(x_size);
 end
 
 Theta = reshape( theta (x_size*nb_channel), x_size, nb_channel);
 
 F = opFFT(x_size);

op = @(x,mode) opRC_measur_matrix_internl(x_size, nb_channel, nb_meas, Theta, W_vector, x, mode);


end

function x = opRC_measur_matrix_internl(x_size, nb_channel, nb_meas, Theta, W_vector, x, mode)
    
if (mode == 1)
    
    
    x = reshape (x, x_size, nb_channel);

    x = ifft(W_vector.*fft(x));
    x = sqrt(nb_meas/x_size)*reshape(Theta.*x,x_size/nb_meas,nb_meas, nb_channel);
    x =  real(sum(x,1));
    x = x(:);
   
elseif (mode == 2)
    
    
    x = reshape(x, nb_meas, nb_channel);
    x_temp(1,:,:) = x;
    clear x;
    
    x_temp = (sqrt(nb_meas/x_size))*repmat(x_temp, x_size/nb_meas, 1);
    x_temp = reshape(x_temp, x_size, nb_channel);
    
    x_temp = real(ifft(conj(W_vector).*(fft(Theta.*x_temp))));
    x = x_temp(:);
  
  
elseif (mode == 0)
     x = {nb_meas,x_size,[0,1,0,1],{'RC_measur_matrix',nb_meas}};
     
else
    fprintf('There is no mode associated to the requested one')
end

end
    