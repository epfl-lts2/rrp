function srec = ai_audio_inpaint_multiple_holes(shole, fs, hole_interval ,param)
%AI_AUDIO_INPAINT_MULTIPLE_HOLES In-paint an audio file with multiples holes
%   Usage:  srec = ai_audio_inpaint_multiple_holes(shole,hole_interval,fs,param)
%   
%   Input parameters:
%       shole       : signal to be in-painted
%       fs          : sampling frequency
%       hole_interval: hole intervals (in seconds, k x 2 matrix) 
%       param       : structure of optional parameters
%   Output parameters:
%       srec        : recomposed signal
%  
%   This function inpaint multiples holes in an audio signal.
%

% Authors: Nathanael Perraudin, Nicki Hollighaus
% Date   : June 2016
 
if nargin<4
    param = ai_conf();
end
 
srec = shole;
hole_interval_fin = hole_interval;
ns_old = size(srec,1);
for ii = 1:size(hole_interval,1)
%     sthole = hole_interval_fin(ii,1); % start of the hole
%     finhole = hole_interval_fin(ii,2); % final of the hole
    param.exclude = hole_interval_fin((ii+1):end,:);
    srecs = ai_audio_inpaint(srec, fs, hole_interval_fin(ii,:), param);
    srec = srecs{1};
 
    %% Adjust hole positions
    [hole_interval_fin, ns_old ] = ai_adjust_hole_position(hole_interval_fin,size(srec,1),ns_old,hole_interval_fin(ii,2),fs);
    
    if param.verbose
        %% Plot some stuff
        plot_tmp(shole,srec,fs, hole_interval, hole_interval_fin)
    end
end
 
 
end
 
function [hole_interval, ns_old ] = ai_adjust_hole_position(hole_interval,ns,ns_old,finhole,fs)
 
for ii = 1:size(hole_interval,1)
    if hole_interval(ii,1)> finhole
        hole_interval(ii,:) = hole_interval(ii,:) + (ns - ns_old)/fs;
    end
end
ns_old = ns;
 
end
 
 
function plot_tmp(shole,srecs,fs, hole_interval, hole_interval_fin)
 
    time0 = (1:size(shole,1)*1.5)/fs;
    
    figure(100)
    subplot(211)
    hold off
    sholet = zeros(numel(time0),size(shole,2));
    sholet(1:numel(shole)) = shole;
    plot(time0,sholet);
    axis tight
    % Plot intervals
    hold on
    plot(time0,interval2signal(hole_interval,time0),'k-');
    
    subplot(212)
    hold off
    srecst = zeros(numel(time0),size(srecs,2));
    srecst(1:numel(srecs)) = srecs; 
    plot(time0,srecst);
    axis tight
    hold on
    plot(time0,interval2signal(hole_interval_fin,time0),'k-');
    
    drawnow
 
end
 
 
function [s ] = interval2signal(hole_interval,time0)
 
nel = size(hole_interval,1);
s = zeros(nel,numel(time0));
 
for ii = 1:nel
    s(ii,:) = double(time0 > hole_interval(ii,1) & time0 < hole_interval(ii,2));
end
 
end

