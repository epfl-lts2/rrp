function [ s ] = check_sound( s )
%CHECK_SOUND check if the sound-signal is mono
%   Usage:  s = check_sound( s );
%
%   Check if the sound is mono (one channel). If not (stereo or
%   multichanel), this function will try to convert it by a simple average.

n_channel_sound=size(s,2);
if n_channel_sound>1
   fprintf('\n   * Warning, sound is not mono. Try to convert it:  ');
   s=mean(s,2);
   fprintf('OK.\n');
end


end

