%AI_DEMO_INPAINT Demonstration of audio inpainting
%
%   Run this script, select your audio file and follow the instruction in
%   the command prompt. Your inpainted file will be saved into the folder
%   'results'.
%
%   In order to run this script you need to install the LTFAT and the
%   GSPBox toolboxes. You can download them at:
%
%   * http://ltfat.github.io/
%   * https://lts2.epfl.ch/gsp/
%
%   In order to gain execution speed, please consider using FLANN with the
%   GSPBox.
%
%   This script is the implementation used in the following contribution:
%   * Paper   : Audio inpainting with similarity graphs
%   * Authors : Nathanael Perraudin, Nicki Holighaus, Piotr Majdak, Peter Balazs 
%   * Date    : June 2016 
%   * ArXiv   : http://arxiv.org/abs/1607.06667
%   * Online demo : https://lts2.epfl.ch/web-audio-inpainting/
%
%   Abstract
%   --------
%   
%   In this contribution, we present a method to compensate for long
%   duration data gaps in audio signals, in particular music. To achieve
%   this task, a similarity graph is constructed, based on a short-time
%   Fourier analysis of reliable signal segments, e.g. the uncorrupted
%   remainder of the music piece, and the temporal regions adjacent to the
%   unreliable section of the signal. A suitable candidate segment is then
%   selected through an optimization scheme and smoothly inserted into the
%   gap.
%
%   References: perraudin2016audio
%


% Authors: Nathanael Perraudin, Nicki Hollighaus
% Date   : June 2016

clear;
close all;
ai_start();


%% Problem creation

[filename,pathname] = uigetfile('*.wav;*.mp3;*.flac','Choose audio file...');   
signame = fullfile(pathname,filename);

[s,fs] = ai_loadsignal(signame);

% hole_interval = [60,62]; % Location of the hole in seconds

while 1
    holeStart = input(sprintf('Enter hole start in seconds (signal length is %.2f s): ',numel(s)/fs));
    if ~isempty(holeStart)
        hole_interval(1) = holeStart;
        break;
    else
        disp('Please try again');
    end
end
while 1
    holeEnd = input(sprintf('Enter hole end in seconds (signal length is %.2f s): ',numel(s)/fs));
    if ~isempty(holeEnd)
        hole_interval(2) = holeEnd;
        if hole_interval(2)<=hole_interval(1)
            fprintf('Please enter position after %.2f s\n',hole_interval(1));
            continue;
        end
        break;
    else
        disp('Please try again');
    end
end

% Create the signal with the hole
shole = s;
shole(hole_interval(1)*fs:hole_interval(2)*fs) = 0;

%% Parameters
% Look inside the ai_conf to change parameters.
param = ai_conf(); 

param.verbose = 1; % Increase the verbosity.

%%

[srecs, G, transitions, info] = ai_audio_inpaint(shole, fs, hole_interval, param);
param = info.param;
%%
if param.verbose
    fprintf('-- Timings of the algorithm --\n')
    fprintf('   a) Features extraction   : %.3g seconds\n',info.timing(1));
    fprintf('   b) Graph construction    : %.3g seconds\n',info.timing(2));
    fprintf('   c) Transitions selection : %.3g seconds\n',info.timing(3));
    fprintf('   d) Signal reconstruction : %.3g seconds\n',info.timing(4));
    fprintf('                              -----\n');
    fprintf(' Total                      : %.3g seconds\n',info.timing(5));
end

%% Save the result and plot the graph
sname = ['experiment','_',num2str(hole_interval(1)),'_',num2str(hole_interval(2))];
srec = srecs{1};

audiowrite2(srec,fs,...
    ['results/',sname,'_rec.wav'])


audiowrite2(shole,fs, ...
    ['results/',sname,'_hole.wav'])

audiowrite2(s,fs,...
    ['results/',sname,'_ori.wav'])
%%
figure(1)
ai_plot_graph(G,transitions{1},hole_interval,param)



