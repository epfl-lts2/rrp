function [D, SR] = audioread2(FN)

% TARGETSR = 0;
% FORCEMONO = 0;
% START = 0; 
% DUR = 0; 
% END = START+DUR;



if exist(FN, 'file') == 0
  error(['audioread2: file ',FN,' not found']);
end

[pth,nam,ext] = fileparts(FN);
ext = lower(ext);

ismp3 = strcmp(ext,'.mp3');
ism4a = strcmp(ext, '.m4a') || strcmp(ext, '.aac') || strcmp(ext, '.mp4');
iswav = strcmp(ext, '.wav');
isrds = strcmp(ext, '.wv1') || strcmp(ext, '.wv2') ...
        || strcmp(ext, '.wvs') || strcmp(ext, '.shn') ...
        || strcmp(ext, '.aif') || strcmp(ext, '.aiff') ...
        || strcmp(ext, '.sph');
isflac = strcmp(ext, '.flac');


if ismp3 || ism4a || iswav || isrds || isflac
  if ismp3
    [D, SR] = mp3read(FN);
  elseif ism4a
    [D, SR]= m4aread(FN);
  elseif iswav
    [D, SR] = wavload(FN);
  elseif isflac
    [D, SR] = flacread(FN);
  else
    error('Unsupported filetype');
  end
else
	error(['audioread2: cannot figure type of file ',FN]);
end
