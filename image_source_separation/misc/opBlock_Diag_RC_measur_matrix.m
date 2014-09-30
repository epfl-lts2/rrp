function op = opBlock_Diag_RC_measur_matrix(x_size, nb_channel, nb_meas, s_gen, t_gen)

%
% OPDLOCK_DIAG_RC_MEASUR_MATRIX creates an operator of Distict Block Diagonal Random convolution measurement
% matrix
%
%   x_size:         input size
%   nb_channel:     number of channels
%   nb_meas:        number of measures
%
% Paper:
%           Compressive sensing 
% Author: Mahdad Hosseini Kamal, Mohammad Golbabaee; ,  Perraudin Nathana?l
% Contact: mahdad.hosseinikamal@epfl.ch, mohammad.golbabaei@epfl.ch,
% nathanael.perraudin@epfl.ch
% Date: June 2012


 W_vector = repmat(s_gen,1,nb_channel);
 
 
 Theta = repmat( t_gen, 1, nb_channel);


op = @(x,mode) opRC_measur_matrix_internl(x_size, nb_channel, nb_meas, Theta, W_vector, x, mode);


end

function x = opRC_measur_matrix_internl(x_size, nb_channel, nb_meas, Theta, W_vector, x, mode)
    
if (mode == 1)
    
    
    x = reshape (x, x_size, nb_channel);

    x = OP_IFFT2_p(W_vector.* OP_FFT2_p(x));
    x = sqrt(nb_meas/x_size)*reshape(Theta.*x,x_size/nb_meas,nb_meas, nb_channel);
    x =  real(sum(x,1));
    x = x(:);
    

elseif (mode == 2)
    
    
    x = reshape(x, nb_meas, nb_channel);
    x_temp(1,:,:) = x;
    clear x;
    
    x_temp = (sqrt(nb_meas/x_size))*repmat(x_temp, x_size/nb_meas, 1);
    x_temp = reshape(x_temp, x_size, nb_channel);
    
    x_temp = real(OP_IFFT2_p(conj(W_vector).*(OP_FFT2_p(Theta.*x_temp))));
    x = x_temp(:);
  
  
elseif (mode == 0)
     x = {nb_meas,x_size,[0,1,0,1],{'RC_measur_matrix',nb_meas}};
     
else
    fprintf('There is no mode associated to the requested one')
end

end