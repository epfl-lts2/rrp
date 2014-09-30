%RR_ALPHA Sensitivity of the alpha parameter
%
%   Reproducible research addendum for phase recovery problem
%   ---------------------------------------------------------
%   
%   AN EXTENDED GRIFFIN LIM ALGORTITHM
%
%   Paper: Perraudin Nathanael, Balazs Peter
%   
%   Demonstration matlab file:  Perraudin Nathanael
%
%   ARI -- April 2013
%   
%   Dependencies
%   ------------
%
%   In order to use this matlab file you need the LTFAT toolbox. You
%   can download it on http://ltfat.sourceforge.net
%   
%   The problem
%   -----------
%
%   From a spectrogram *S*, we would like to recover the signal with the
%   closest spectrogram.
%
%   The problem could be written in the form:
%
%   ..   minimize_x   || |Gx| - S ||_2
%
%   .. math:: \| |Gx| - S \|_2
%
%   with 
%         - $S$ :  The original spectrogram
%         - $G$ :  A frame
%         - $x$ :  The signal we would like to recover
%
%   Note that for these simulations, we know that it exist one x such that
%   $|Gx|=S$. 
%
%
%   Algorithms for solving the problem
%   ----------------------------------
%
%   We will compare 2 different algorithms to solve the problem  
%
%   * Griffin-Lim : Original algorithm designed by griffin and Lim
%
%   * Forward-PBL : A modification of the Griffin-Lim
%
%   
%   Goal of these simulations
%   -------------------------
%
%   Observe the role of the parameter 'alpha' in th Foward-PBL algorithm.
%   
%   * For 'alpha'=0, we recover the Griffin lim algorithm
%   * For 'alpha'>1, the algorithm is usually unstable
%   * The optimum seems to be 'alpha'close to 1.
%   
%
%   Results
%   -------
%
%   .. figure::
%
%      Different values of 'alpha' on 'gspi'
%
%       
%
%   .. figure::
%
%      Different values of 'alpha' on 'bat'
%
%       
%


% ----------------------------------------------------------------------- %
% Author: Nathanael Perraudin
% April 2013
% FLAME project, ARI (Acoustic research institute), Vienna
% ----------------------------------------------------------------------- %

%% Initialisation
clear all;
close all;
clc;

% The LTFAT toolbox is required for this demonstration file
ltfatstart % if this line is a problem, add to path the LTFAT toolbox.

%% General Parameters
alpha=[0.2 0.7 0.9 0.99 1 1.05 1.2]; % Parameter for the proposed algorithm
                                
% Frame parameter
a=32;                               % Size of the shift in time
M=256;                              % Number of frequencies     
window_type='hann';                % Type of window

real_signal=1;                      % we work with real signal only
type_multiplier='full';


% algorithm parameter

solving_method={'GLA','FGLA'};   

param.verbose=0;                    % Display parameter

% Starting point
param.zero_phase=1;                 % Starting with a zero phase
param.random_phase=0;               % Starting with a random phase

param.tol=0;                        % Tolerance to stop iterating 

% Plot parameter
paramplot.pathfigure='comparaison/';
paramplot.legendlocation='NorthWest';
paramplot.save = 0;
paramplot.position=[100 100 600 400];

%% Lauch different simulations with different signals

%% gspi
param.maxit=10000;                  % Maximal number of iteration
sound_name='bat';                  % Sound name
main(a,M,alpha,window_type,sound_name,real_signal,type_multiplier,...
                                      solving_method,param,paramplot);

                                      
%% bat       
param.maxit=1000;                  % Maximal number of iteration
sound_name='gspi';                   % Sound name                   
main(a,M,alpha,window_type,sound_name,real_signal,type_multiplier,...
                                      solving_method,param,paramplot); 
                                  