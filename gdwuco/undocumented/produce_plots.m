% Produce all plots for the article

global GLOBAL_save
global GLOBAL_baw
global GLOBAL_figpaper

GLOBAL_baw = 0;         % to save the result in black and white
GLOBAL_save = 1;        % or 1 if you want to save the figures
GLOBAL_figpaper= 1;     % to make the figure for the paper 

rr_introduction;
rr_beat_itersine;
rr_gabfirdual_smoothness;
rr_gabfirtight;
rr_gaboptdual;
rr_tradeoff_smoothness;
rr_gabfirdual

%%
GLOBAL_figpaper= 0;     % to make the figure for the paper 

rr_gabfirtight;
rr_beat_itersine;
