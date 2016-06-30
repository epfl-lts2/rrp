function [  ] = ai_start(  )
%AI_START Add the necessary path for the audio in-painting
%   Usage: ai_start();
%   

ltfatstart;
gsp_start;

path = fileparts(mfilename('fullpath'));
addpath(genpath(path));

end

