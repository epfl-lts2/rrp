%RR_SIGNALS Phase recovery problem onto different signals
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
%   Test the algorithm on various acoustic signals.
%
%   We use all the test signal from the toolbox LTFAT. For each of them the
%   Foward-PBL is better than the Griffin-Lim.
%
%
%   Results
%   -------
%
%   .. figure::
%
%      Phase recovery problem on 'gspi'
%
%       
%
%   .. figure::
%
%      Phase recovery problem on 'traindoppler'
%
%       
%
%   .. figure::
%
%      Phase recovery problem on 'cocktailparty'
%
%       
%
%   .. figure::
%
%      Phase recovery problem on 'linus'
%
%       
%
%   .. figure::
%
%      Phase recovery problem on 'bat'
%
%       
%
%   .. figure::
%
%      Phase recovery problem on 'greasy'
%


% ----------------------------------------------------------------------- %
% Author: Nathanael Perraudin
% April 2013
% FLAME project, ARI (Acoustic research institute), Vienna
% ----------------------------------------------------------------------- %

%% Initialisation
clear;
close all;

% The LTFAT toolbox is required for this demonstration file
ltfatstart % if this line is a problem, add to path the LTFAT toolbox.

%% General Parameters
                                
% Frame parameter
a=32;                               % Size of the shift in time
M=256;                              % Number of frequencies     
window_type='nuttall';                % Type of window

real_signal=1;                      % we work with real signal only
type_multiplier='full';


% algorithm parameter

solving_method={'GLA','FGLA'};   

param.maxit=1000;                  % Maximal number of iteration

param.verbose=0;                    % Display parameter

% Starting point
param.zero_phase=1;                 % Starting with a zero phase
param.random_phase=0;               % Starting with a random phase

param.tol=0;                        % Tolerance to stop iterating 

alpha=0.99;                         % Parameter for the proposed algorithm

% Path for saving the figures
paramplot.pathfigure='comparaison/';
paramplot.legendlocation='SouthEast';
paramplot.save = 0;
paramplot.position=[100 100 300 225];

%% Lauch different simulations with different signals
                                        

%% gspi
sound_name='gspi';                  % Sound name
main_flga(a,M,alpha,window_type,sound_name,real_signal,type_multiplier,...
                                      solving_method,param,paramplot);

                                        
%% traindoppler                                        
sound_name='traindoppler';          % Sound name
main_flga(a,M,alpha,window_type,sound_name,real_signal,type_multiplier,...
                                      solving_method,param,paramplot);  
                                        
%% cocktailparty                                        
sound_name='cocktailparty';          % Sound name
main_flga(a,M,alpha,window_type,sound_name,real_signal,type_multiplier,...
                                      solving_method,param,paramplot);
                                        
%% linus                                        
sound_name='linus';                 % Sound name
main_flga(a,M,alpha,window_type,sound_name,real_signal,type_multiplier,...
                                      solving_method,param,paramplot);
                                        
%% bat                                        
sound_name='bat';                   % Sound name
param.maxit=20000;                    
main_flga(a,M,alpha,window_type,sound_name,real_signal,type_multiplier,...
                                      solving_method,param,paramplot);                                         
                                        
%% greasy     
param.maxit=10000;
sound_name='greasy';                % Sound name
main_flga(a,M,alpha,window_type,sound_name,real_signal,type_multiplier,...
                                      solving_method,param,paramplot);                                       
