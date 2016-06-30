function [featuremat, param, patchmax] = ai_compute_features(s,fs,param)
%AI_COMPUTE_FEATURES compute a time-frequency features matrix
%   Usage:  featuremat = ai_compute_features(s,fs)
%           featuremat = ai_compute_features(s,fs,param)
%           [featuremat, param] = ai_compute_features(...)
%           [featuremat, param, patchmax] = ai_compute_features(...)
%
%   Input parameters:
%       s           : audio signal (vector)
%       fs          : sampling frequency
%       param       : structure of optional parameters
%   Output parameters:
%       featuremat  : matrix of features
%       param       : updated parameters
%       patchmax    : vector of maximum values of the Gabor transform
%
%   *param.featuretype* can be set to an integer to change the type of
%   features computed.
%
%   * 0: log spectrogram
%   * 1: spectrogram
%   * 2: reassigned spectrogram
%   * 3: phase derivatives
%   * 4: spectrogram + phase derivatives
%   * 5: reassigned specectrogram + phase derivatives
%   * 6: MFCC
%
%   *param.verbose* can be set to *1* to plot the feature matrix
%
%   This function requires the LTFAT toolbox to work.
%
 
 
% Authors: Nicki Hollighaus, Nathanael Perraudin
% Date   : June 2016
 
if nargin<3
    param = struct;
end
 
if ~isfield(param, 'featuretype')
    param.featuretype = 4;
end
 
if ~isfield(param, 'verbose')
    param.verbose = 0;
end
if ~isfield(param,'win_length'), param.win_length = 1024; end
if ~isfield(param,'win'), param.win = 'itersine'; end
if ~isfield(param,'M'), param.M = param.win_length ; end
if ~isfield(param,'a'), param.a = 128 ; end
if ~isfield(param,'lambda'), param.lambda = 3/2; end
if ~isfield(param,'cepstralcoefs'), param.cepstralcoefs = 25 ; end
if ~isfield(param,'dbrange'), param.dbrange = 50 ; end
 
 
% TODO: I think this is not good!
% Make the signal a suitable length
Loriginal_r = length(s);
L=dgtlength(param.win_length,param.a,param.M);
s = s( 1: (Loriginal_r - mod(Loriginal_r,L)));
 
g = firwin(param.win,param.win_length);
Gd = @(x) dgtreal(x,g,param.a,param.M);
 
F = frame('dgt',g,param.a,param.M);
[~,B] = framebounds(F);
 
 
gabor_transform = Gd(s)/sqrt(B);
mxGT = max(abs(gabor_transform(:)));
[vsize, owidth]= size(gabor_transform);     
%vsize = vsize/2+1;
 
% Compute feature matrix
 
switch param.featuretype
    case {4,5}
        num_features = 2;
    otherwise
        num_features = 1;
end
 
if param.featuretype == 6
    vsizefull = param.cepstralcoefs+(num_features-1)*vsize-1;
else
    vsizefull = num_features*vsize;
end
 
featuremat = zeros(vsizefull,owidth);
 
switch param.featuretype
    case {2,3,4,5} % Precompute phase derivatives
        [tgrad,fgrad] = gabphasegrad('dgt',s,g,param.a,param.M);
        tgrad = tgrad(1:floor(end/2)+1,:);
        fgrad = fgrad(1:floor(end/2)+1,:);
    otherwise        
end
 
patchmax = ones(1,owidth);
if nargout > 1
    patchmax = zeros(1,owidth);
    for kk = 1:owidth
       patchmax(kk) = norm(gabor_transform(:,kk),Inf);       
    end
    patchmax(patchmax < 10e-8) = Inf;
end
 
switch param.featuretype
    case {0,4} %log spectrogram feature
        temp = abs(gabor_transform);
        temp = 20*log10(temp);
        temp = temp - max(temp(:));
        temp = temp + param.dbrange;
        temp(temp<0) = 0;
        featuremat(1:vsize,:) = temp/max(temp(:));
    case {1} % spectrogram feature
        temp = abs(gabor_transform).^2;
        temp = temp./repmat(patchmax.^2,vsize,1);
        featuremat(1:vsize,:) = temp;
    case {2,5} % reassigned spectrogram feature        
        temp = gabreassign(abs(gabor_transform).^2,tgrad,fgrad,param.a);
        temp = 10*log10(temp);
        temp = temp - max(temp(:));
        temp = temp + 1.5*param.dbrange;
        temp(temp<0) = 0; 
        %temp = temp./repmat(patchmax.^2,vsize,1);
        featuremat(1:vsize,:) = temp/max(temp(:));
    case {6}   % MFCC feature
%        error('This will be implemented later on');
        temp = mfcc_from_dgtreal(gabor_transform,fs,param.a,param.M,param.cepstralcoefs);
        temp = abs(temp(2:end,:));
        temp = temp/max(temp(:));
        temp = temp./repmat(patchmax.^2,param.cepstralcoefs,1);
        featuremat(1:param.cepstralcoefs-1,:) = temp;
    otherwise        
end
 
switch param.featuretype
    case {3,4,5} % phase derivative feature
       temp = tgrad(1:vsize,:);
       temp1 = (abs(gabor_transform(1:vsize,:))> 2*mxGT*10e-3);
       temp = temp1.*temp;
       
       clear temp1; 
       
       kernel = fir2long(firwin('hann',8),size(temp,2));
       temp = pconv(temp.',kernel).';
%        for kk = 1:size(temp,1)
%           temp(kk,:) = pconv(temp(kk,:),kernel);
%        end
       mx = max(abs(temp(:)));
       temp = temp/mx/param.lambda;
       featuremat((num_features-1)*vsize+(1:vsize),:) = temp;
    otherwise        
end
 
% else
%     % Old method
%     
% 
%     
%     g = firwin(param.win,param.win_length);
%     Gd = @(x) dgtreal(x,g,param.a,param.M);
%     F = frame('dgt',g,param.a,param.M);
%     [~,B] = framebounds(F);
%     gabor_transform = Gd(s)/sqrt(B);
%     [vsize, owidth]= size(gabor_transform);     
% 
%      patchmax = ones(1,owidth);
%     if nargout > 1
%         patchmax = zeros(1,owidth);
%         for kk = 1:owidth
%            patchmax(kk) = norm(gabor_transform(:,kk),Inf);       
%         end
%         patchmax(patchmax < 10e-8) = Inf;
%     end
%     spec = abs(gabor_transform).^2;
%     transpec = 20*log10(spec);
%     transpec = transpec - max(transpec(:));
%     transpec = transpec + param.dbrange;
%     transpec(transpec<0) = 0;
%     featuremat = transpec/max(transpec(:));
% end
 
if param.verbose>1
    if param.featuretype == 4;
        figure(101)
        subplot(211)
        time = (1:round(vsize/4))/round(vsize/4)*length(s)/fs;
        frequency = (0:size(featuremat,2)-1)/size(featuremat,2)*fs/2;
        imagesc(time,frequency,abs(featuremat(1:round(vsize/4),:)))
        axis('xy');
        xlabel('Time (s)')
        ylabel('Frequency (Hz)')
        colorbar;
        title('Normalized log-spectrogram')        
        subplot(212)
        imagesc(time,frequency,featuremat(vsize+1:vsize+round(vsize/4),:))
        axis('xy');
        xlabel('Time (s)')
        ylabel('Frequency (Hz)')
        colorbar;
        title('Phase derivative')   
    else
        figure(101)
       imagesc(abs(featuremat))
        axis('xy');
        xlabel('Patch number')
        ylabel('Frequency feature number')
        colorbar;
        title('Feature space')
    end
end
 
if nargout > 1
    % Why to we do this?
    patchmax(patchmax > 10e10) = 0;
end
 
end

