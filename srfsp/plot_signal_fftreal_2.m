function [ s ] = plot_signal_fftreal_2( x,p,fig,fs,sol )
%PLOT_SIGNAL_FFTREAL Plot image plugin for the UNLocBoX
%   Usage [ s ] = plot_signal_fftreal( im,fig );
%
%   Input parameters:
%         x     : Structure of data
%         p     : Number of iteration between 2 plots...
%         fig   : Figure
%         fs    : Sampling frequency
%
%   Output parameters:
%         s     : Input image
%
%   This plugin display the signal every iterations of an algorithm. To use
%   the plugin juste define::
%       
%       fig=figure(100);
%       param.do_sol=@(x) plot_signal_fftreal(x,p,fig);
%
%   In the structure of optional argument of the solver.
%

% Author: Nathanael Perraudin
% Date  : 3rd april 2014

if ~mod(x.iter-1,p)
    % select the figure
    if x.iter<2
        figure(fig);
    end
    % display the signal
    fftx = fftreal(x.sol(:));
    N = length(fftx);

    w = linspace(0,fs-fs/N,N);
    subplot(121)
    plot(w,abs(fftx));
    xlabel('Frequency')
    ylabel('Modulus')
    title(['Current it: ', num2str(x.iter),'   Curr obj: ', ...
        num2str(x.info.objective(x.iter))]);
    subplot(122)
    plot(w,abs(sol(:)));
    xlabel('Frequency')
    ylabel('Modulus')
    title(['Solution,  diff = ',num2str(norm(fftx(:)-sol(:)))]);
    drawnow;
end


s=x.sol;



end
