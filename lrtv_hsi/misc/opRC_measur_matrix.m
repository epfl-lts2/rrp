function op = opRC_measur_matrix(x_size,nb_meas)

%
% OPRC_MEASUR_MATRIX creates an operator of Random convolution measurement
% matrix
%
%   x_size:         input size
%   nb_meas:        number of measures
%
% Paper:
%           Compressive sensing 
% Author: Mahdad Hosseini Kamal
% Contact: mahdad.hosseinikamal@epfl.ch
% Date: Oct. 2010


 W_vector = sigma_gen(x_size);
 Theta = theta (x_size);
 


op = @(x,mode) opRC_measur_matrix_internl(x_size, nb_meas, Theta, W_vector, x, mode);


end

function y = opRC_measur_matrix_internl(x_size, nb_meas, Theta, W_vector, x, mode)
    
if (mode == 1)


    H = ifft(W_vector.*fft(x));
    Bk = sqrt(nb_meas/x_size)*reshape(Theta.*H,x_size/nb_meas,nb_meas);
    y_temp =  real(sum(Bk,1));
    y = y_temp;
    y = y(:);
   
elseif (mode == 2)
    

    
    x_temp = (sqrt(nb_meas/x_size))*repmat(transpose(x(:)), x_size/nb_meas, 1);
    x_temp = x_temp(:);
    
    y_temp = real(ifft(conj(W_vector).*(fft(Theta.*x_temp))));
    y = y_temp;
    y = y(:);
  
  
elseif (mode == 0)
     y = {nb_meas,x_size,[0,1,0,1],{'RC_measur_matrix',nb_meas}};
     
else
    fprintf('There is no mode associated to the requested one')
end

end
    