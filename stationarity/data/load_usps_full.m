function [x, y, xx, yy] = load_usps_full()
%LOAD_USPS_FULL load the USPS dataset
%   Usage: [x, y, xx, yy] = load_usps_full();
%
%   Ouptut parameters
%       x       : training images
%       y       : training labels
%   Output parameters
%       xx      : testing images
%       yy      : testing labels


try
  load usps_resampled.mat
catch
  error('Error: the file usps_resampled.mat was not found.')
end

y = matrix2label(train_labels',0);
x = train_patterns;

yy = matrix2label(test_labels',0);
xx =  test_patterns;

end