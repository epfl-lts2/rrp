function [ s,fs ] = load_sound( name )
%LOAD_SOUND Load one of the following sound: 
%   Usage: s = load_sound( name );
%          [ s,fs ] = load_sound( name );
%
%   Input parameters:
%         name  : Name of the sound.
%   
%   Output parameters:
%         s     : signal.
%         fs    : sampling frequency.
%
%   This function load a sound from the LTFAT toolbox by is name. The name
%   availlable are:
%
%   * bat
%   * cocktailparty
%   * greasy
%   * gspi
%   * linus
%   * otoclick
%   * traindoppler
%
%   In order to speed up computation it return only "power of 2" numbers of
%   samples.
%

% Load the sound
if strcmp(name,'bat')
   [s,fs]= bat();
elseif strcmp(name,'cocktailparty')
   [s,fs]= cocktailparty();
elseif strcmp(name,'greasy')
   [s,fs]= greasy();
elseif strcmp(name,'gspi')
   [s,fs]= gspi();
elseif strcmp(name,'linus')
   [s,fs]= linus();
elseif strcmp(name,'otoclick')
   [s,fs]= otoclick();
elseif strcmp(name,'traindoppler')
   [s,fs]= traindoppler();
else
   error('Incorrect sound name!'); 
end



% Cut the extra part of th sound. We work with "power of 2" number of
% samples.
%length(s);
%N=floor(log(length(s))/log(2));
%s=s(1:2^N);

end

