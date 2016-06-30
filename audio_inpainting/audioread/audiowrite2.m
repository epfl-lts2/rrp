function audiowrite2(D,SR,FN)

try
    audiowrite(FN,D,SR);
catch

    [pth,nam,ext] = fileparts(FN);

    % Ensure the output directory exists
    if ~exist(pth, 'dir')
      mkdir(pth);
    end

    ext = lower(ext);
    if strcmp(ext,'.mp3')
      mp3write(D,SR,FN);
    elseif strcmp(ext, '.m4a') || strcmp(ext, '.aac') || strcmp(ext, '.mp4')
      m4awrite(D,SR,FN);
    elseif strcmp(ext, '.flac')
      flacwrite(D,SR,FN);
    elseif strcmp(ext, '.wv1') || strcmp(ext, '.wv2') || strcmp(ext, '.shn')
      % nist sphere
      % use dpwelib as a mex file
      audiowrite(D,SR,FN);
    else
      %wavwrite(D,SR,FN);
      wavsave(D,SR,FN);
    end

end
