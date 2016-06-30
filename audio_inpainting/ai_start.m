function [  ] = ai_start(  )
%AI_START Start the audio in-painting box
%   Usage: ai_start();
%   
%   This function add the necessary path for the audio in-painting.
%

ltfatstart;
gsp_start;

path = fileparts(mfilename('fullpath'));
addpath(genpath(path));

end

