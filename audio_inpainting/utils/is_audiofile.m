function [bool, name, pth] = is_audiofile(FN)

[pth,nam,ext] = fileparts(FN);
ext = lower(ext);

ismp3 = strcmp(ext,'.mp3');
%ism4a = strcmp(ext, '.m4a') || strcmp(ext, '.aac') || strcmp(ext, '.mp4');
iswav = strcmp(ext, '.wav');
% isrds = strcmp(ext, '.wv1') || strcmp(ext, '.wv2') ...
%         || strcmp(ext, '.wvs') || strcmp(ext, '.shn') ...
%         || strcmp(ext, '.aif') || strcmp(ext, '.aiff') ...
%         || strcmp(ext, '.sph');
isflac = strcmp(ext, '.flac');
name = [pth, '/',nam];

bool = ismp3 + iswav + isflac;

end

