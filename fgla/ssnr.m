function [ ssnrm ] = ssnr( S_ref,S_est )
%SSNR Spectrogram signal to noise ratio
%   Usage: [ ssnrm ] = ssnr( S_ref,S_est )
%
%   Input parameters:
%         S_ref : reference spectrogram.
%         S_est : estimate spectrogram.
%   
%   Output parameters:
%         ssnrm : ssnr measurement.


ssnrm=snr(abs(S_ref),abs(S_est));
%fsnr=sum(sum(abs(log(abs(Sref)+eps)-log(abs(Sest)+eps))));

end

function snr = snr(s_in, s_est)
%SNR Compute the SNR between two signals in dB
%   Usage X = snr(s_in, s_est)
 
snr = norm(s_in(:)-s_est(:))/norm(s_in(:));
snr = - 10 * log10(snr);

end
