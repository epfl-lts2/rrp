function [ s ] = plot_window( x,p,fig,filename)
%PLOT_WINDOW Plot image plugin for the UNLocBoX
%   Usage [ s ] = plot_window( x,p,fig,filename);
%             s = plot_window( x,p,fig);
%             s = plot_window( x,p);
%             s = plot_window( x);
%
%   Input parameters:
%         x     : Structure of data
%         p     : Number of iteration between 2 plots...
%         fig   : Figure
%         filename: filename for the gif
%
%   Output parameters:
%         s     : Input image
%
%   This plugin display the window every iterations of an algorithm. To use
%   the plugin juste define::
%       
%       fig=figure(100);
%       param.do_sol=@(x) plot_window(s,p,fig);
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
    g = real(x.sol);
    Lg = length(g);


    if mod(Lg,2)
        a =-(Lg-1)/2:(Lg-1)/2;
    else
        a = -Lg/2:(Lg/2-1);
    end

    plot(a,fftshift(g));
    xlim([min(a),max(a)]);
    title(['Current it: ', num2str(x.iter),'   Curr obj: ', ...
        num2str(x.info.objective(x.iter))]);

    drawnow;
end

% return the signal
s=x.sol;

if filename
      frame = getframe(fig);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);  
      dt = 0.1;
    if x.iter == 2;
        imwrite(imind,cm,[filename,'.gif'],'gif','DelayTime',dt,...
             'Loopcount',inf,'writemode','overwrite');
    else
    	imwrite(imind,cm,[filename,'.gif'],'gif','DelayTime',dt,...
            'WriteMode','append');
    end
end



end

