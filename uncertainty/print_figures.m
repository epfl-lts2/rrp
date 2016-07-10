%PRINT_FIGURES print all the figures for the article
%
%   This script print all the figures for the article ...

% Author: Nathanaël Perraudin
% Date  : 

clear all; %#ok<CLALL>
close all;

global SAVE ;
global sf;
global type;


% set the variable to save the results
SAVE = 1;

% launch all the scripts
global_illustration;
sparsity_illustration;
comet_coherence;
modified_path_eigenvectors;
modified_path_coherence;
norm_of_localization_operator;
uncertainty_bounds;
spread_gabor;
sf = 0; %#ok<NASGU>
local_uncertainty;
sf = 1;
local_uncertainty;
atom_localization;
modified_path_gabor;

type = 'sensor'; %#ok<NASGU>
adaptated_sampling;

type = 'community';
adaptated_sampling;

