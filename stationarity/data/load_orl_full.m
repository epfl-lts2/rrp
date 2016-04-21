function [x, y, xx, yy] = load_orl_full()
%LOAD_ORL_FULL load the ORL dataset
%   Usage: [x, y, xx, yy] = load_orl_full();
%
%   Ouptut parameters
%       x       : training images
%       y       : training labels
%   Output parameters
%       xx      : testing images
%       yy      : testing labels


try
  load ORL_data.mat
catch
  disp('Error: the file ORL_data.mat was not found. Perhaps you need to')
  disp('run the script get_dataset_...sh')
  x = []; y = []; xx = []; yy = [];
  return
end

x = double(X)';
y = labels(:,1);
xx = [];
yy = [];
end