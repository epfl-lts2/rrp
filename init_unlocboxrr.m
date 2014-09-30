function init_unlocboxrr()
%INIT_TOOLBOXRR Initialize the toolbox
%   Usage: init_unlocboxrr()
%
%   Initialisation script for the convex optimization problems toolbox
%   This script add the different path needed to run the toolbox


% Author: Nathanael Perraudin
% E-mail: nathanael.perraudin@epfl.ch
% Date: jan 2013


%% adding dependency

global GLOBAL_path;
GLOBAL_path = fileparts(mfilename('fullpath'));


addpath(genpath(GLOBAL_path));


