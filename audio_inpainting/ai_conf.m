function param = ai_conf()
% Global parameters for the algorithm
 
% Authors: Nathanael Perraudin, Nicki Hollighaus
% Date   : June 2016
 
% Configuration parameters
param.featuretype = 4;  % 0 - log spectrogram
                        % 1 - spectrogram, 2 - log reassigned spectrogram, 
                        % 3 - phase derivatives, 4 - log spec + phase,
                        % 5 - log reass spec + phase, 6 - MFCC
                        
param.finetune = 'wave'; % 'none' - no fine-tuning                        
                         % 'wave' - waveform correlation
                         % 'beat' - beat matching
                        
param.lambda = 3/2;       % trade-off between the features    
param.cepstralcoefs = 25;
param.a = 128;          % Hope in time
param.dbrange = 50;
param.win = 'itersine';
param.win_length = 1024;
param.M = param.win_length;
 
param.fsmax = 12000;      % Upper bound on the sampling frequency
 
param.k = 40;
param.graph_type = 'knn'; % radius or knn (radius is very slow)
                          % If radius is activated decrease param.k and
                          % increase a bit the threshold
param.threshold = 2;      % Threshold to cut small connection
param.diagdist = 40;      % Size of the convolution kernel
param.use_flann = 1;    % Much faster but a slightly more random
 
param.melt = 2;         % Use melted transitions:
                        % 0 - Do not use, 1 - time domain cross-fading,
                        % 2 - time-frequency cross-fade (this will allow 
                        % for phase adjustments)
 
param.verbose = 1;      % Plot some info figures
 
param.keeplength = 0;   % The recovered signal has the same length as the initial signal
param.keepsignal = 0;   % The algorithm tries to keep as much as possible from the original signal
 
param.ns = 1;           % Number of different reconstructions
 
param.subgraph.time_range = 5; % Subgraph parameter number of seconds before and after the gap.
 
 
 
end

