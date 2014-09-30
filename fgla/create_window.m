function [g,info_g ] = create_window(a,M,L,window_type)
%CREATE_WINDOWS Create a windows
%   Usage: g = create_windows(a,M,L,windows_type)
%          [g,info_g ] = create_windows(a,M,L,windows_type)
%
%   Input parameters:
%         a     : Shift in time.
%         M     : Number of frequency channel.
%         L     : Length of the signal.
%         windows_type: type of window.
%   
%   Output parameters:
%         g     : window.
%         info_g: information about the window.
%
%   This function create a window from the frame parameter. The different
%   type of windows are:
%
%   * gauss         : Gaussian
%   * hann          : Hann
%   * hamming       : Hamming
%   * blackman      : Blackman
%   * square        : Square
%   * tria          : Triangular
%   * nuttall       : Nuttall
%   * itersine      : Itersine (tight)
%   * sqrttria      : Square of triangulare
%

mg=M; % size of the support of the windows. 

if strcmp(window_type,'gauss')
    [g,info_g]=gabwin('gauss',a,M,L);

elseif strcmp(window_type,'hann')
    [g,info_g]=gabwin({'hann',mg},a,M,L);
    
elseif strcmp(window_type,'hamming')
    [g,info_g]=gabwin({'hamming',mg},a,M,L);
    
elseif strcmp(window_type,'blackman')
    [g,info_g]=gabwin({'blackman',mg},a,M,L);
    
elseif strcmp(window_type,'square')
    [g,info_g]=gabwin({'square',mg},a,M,L);
    
elseif strcmp(window_type,'tria')
    [g,info_g]=gabwin({'tria',mg},a,M,L);
    
elseif strcmp(window_type,'nuttall')
    [g,info_g]=gabwin({'nuttall',mg},a,M,L);
    
elseif strcmp(window_type,'itersine')
    [g,info_g]=gabwin({'itersine',mg},a,M,L);
    
elseif strcmp(window_type,'sqrttria')
    [g,info_g]=gabwin({'sqrttria',mg},a,M,L); 
else
   error('Unknown windows type!');
end

    
end

