function [s,fs] =ai_loadsignal(filename,ratio)
%AI_LOADSIGNAL Load a signal, down-sample it, and make it mono
%   Usage:  [s, fs] =ai_loadsignal(filename)
%           [s, fs] =ai_loadsignal(filename, ratio)
%   
%   Input parameters:
%       filename    : name of the file (text)
%       ratio       : down-sampling ratio
%   Output parameters:
%       s           : signal
%       fs          : sampling frequency
%
%   This function load a signal, down-sample it and make it mono.
%


if nargin<2
    ratio = 1;
end

% read the file
try
    [s,fs] = audioread(filename);
catch
    [s,fs] = audioread2(filename);
end

if size(s,2) >=1
    s = mean(s,2);
end
fs = fs/ratio;
s = resample(s,1,ratio);
end