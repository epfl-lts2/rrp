function [srec, G, transitions, info] = ai_audio_inpaint(shole, fs, hole_interval, param)
%AI_AUDIO_INPAINT In-paint an audio file
%   Usage:  srec = ai_audio_inpaint(shole, fs, sthole, finhole);
%           srec = ai_audio_inpaint(shole, fs, sthole, finhole, param);
%           [srec, G] = ai_audio_inpaint(...);
%           [srec, G, jumps] = ai_audio_inpaint(...);
%           [srec, G, jumps, info] = ai_audio_inpaint(...);
%   
%   Input parameters:
%       shole       : signal to be in-painted
%       fs          : sampling frequency
%       hole_inteval : start and end position of the hole (in seconds)
%       param       : structure of optional parameters
%   Output parameters:
%       srec        : recomposed signal
%       G           : graph used for the reconstruction
%       jumps       : positions of the transitions on the graph
%       info        : structure of additional information
%   
%   This function performs in-paint a hole in an audio file. It follows the
%   following steps:
%       
%   1. It creates a graph of transitions of the audio file. See the function
%      |ai_time_audio_graph| for more details.
%   2. It searches for two optimal transitions.
%   3. It reconstructs the signal using these two optimal transitions.
%
%   References: perraudin2016audio
%
 
% Authors: Nathanael Perraudin, Nicki Hollighaus
% Date   : June 2016
 
timetot = tic;
 
if hole_interval(1)<3
    error('The hole cannot be less than 3 seconds after the beginning of the audio file');
end
 
if hole_interval(2) > length(shole)/fs - 3 
    error('The hole cannot be more than 3 seconds before the end of the audio file');
end
 
 
 
sthole = hole_interval(1)*fs; % start of the hole
finhole = hole_interval(2)*fs; % final of the hole
 
if nargin<4
    param = ai_conf();
end
 
if ~isfield(param,'ns'), param.ns = 1; end
if ~isfield(param,'exclude'), param.exclude = [] ; end
if ~isfield(param,'threshold'), param.threshold = 2 ; end
if ~isfield(param,'verbose'), param.verbose = 1 ; end
if ~isfield(param,'keepsignal'), param.keepsignal = 0; end
if ~isfield(param,'keeplength'), param.keeplength = 0; end
if ~isfield(param,'finetune'), param.finetune = 'wave'; end
if ~isfield(param,'melt'), param.melt = 2; end
 
qvecfin = cell(param.ns,1);
 
% Use only the subgraphs
param.subgraph.time_in = sthole;
param.subgraph.time_fin = finhole;
 
 
no_solution_found = 1;
Ntry = 1;
while no_solution_found
    
    %% 1) build the graph
    [G,Gfull,param,timing]  = ai_time_audio_graph(shole,fs,param);
    time = G.time;
    if param.verbose
        disp('Find optimal transitions');
    end
    
    %% 2) find the transitions
    % a) convert from time to node
    timetransition = tic;
    [~,ind1] = min(abs(time-sthole));
    [~,ind2] = min(abs(time-finhole));
 
    % If some part of the signal needs to be excluded
    exclude = zeros(size(param.exclude));
    for ii = 1:size(param.exclude,1)
        for jj = 1:size(param.exclude,2)
            [~,exclude(ii,jj)] = min(abs(time-param.exclude(ii,jj)));
        end
    end
 
    % Search for the optimal transition(s)
    if param.keeplength
        if param.keepsignal
            [transitions, qvec] = ai_find_transitions(Gfull,ind1,ind2,10,100,10,param.ns,exclude );        
        else
            [transitions, qvec] = ai_find_transitions(G,ind1,ind2,1,100,1,param.ns, exclude ); 
        end
    else
        if param.keepsignal
            [transitions, qvec] = ai_find_transitions(Gfull,ind1,ind2,100,100,1,param.ns,exclude);        
        else
            [transitions, qvec] = ai_find_transitions(G,ind1,ind2,1,100,1,param.ns,exclude ); 
        end
    end
    
   if Ntry == 1
        nh = std(shole(:));
        for ii=1:numel(qvec)
            qvecfin{ii} = [qvec{ii}, sqrt(G.sigma), nh, sqrt(G.sigma)/nh, norm(Gfull.W,'fro'),norm(G.W,'fro'),G.Ne ]; %#ok<AGROW>
        end
   end
    
    if isnan(sum(transitions{1}(:)))
        param.threshold = param.threshold/2;
        Ntry = Ntry + 1;
    else
        no_solution_found = 0;
    end
    timing_transition = toc(timetransition);
    
    if Ntry > 5
        error('The algorithm did not find any solutions')
    end
end
 
 
% 3) Reconstruct the signal    
time_reconstruction = tic;
if param.verbose
    disp('Reconstruct the signal');
end
 
%% Compute some numbers (to be updated)
try
    rel_diff = @(x,y) norm(x(:)-y(:))^2 / ( norm(x(:)) * norm(y(:)) );
 
    features = G.features;
    dd = 40;
    convk = tripuls(-dd:dd,2*dd+1); % kernel
    convk1 = convk;
    convk2 = convk;
    convk1((dd+1):end) = 0; 
    convk2(1:(dd-1)) = 0;
 
    mnf = norm(features,'fro')/sqrt(size(features,2))*norm(convk1);
    for ii = 1:length(qvecfin)
        qvecfin{ii}(end+1) = rel_diff( features(:,transitions{ii}(1,1)+(-dd:dd)) * diag(convk1), ...
                        features(:,transitions{ii}(1,2)+(-dd:dd)) * diag(convk1) );
        qvecfin{ii}(end+1) = rel_diff( features(:,transitions{ii}(2,1)+(-dd:dd)) * diag(convk2), ...
                        features(:,transitions{ii}(2,2)+(-dd:dd)) * diag(convk2) );
        qvecfin{ii}(end+1) = ( norm(features(:,transitions{ii}(1,1)+(-dd:dd)) * diag(convk1),'fro') * ...
                      norm(features(:,transitions{ii}(1,2)+(-dd:dd)) * diag(convk1),'fro') / ...
                      mnf^2 )^(-1);
        qvecfin{ii}(end+1) = ( norm(features(:,transitions{ii}(2,1)+(-dd:dd)) * diag(convk2),'fro') * ...
                      norm(features(:,transitions{ii}(2,2)+(-dd:dd)) * diag(convk2),'fro') / ...
                      mnf^2 )^(-1);
    end
end
%%
 
srec = ai_reconstruct_signal(G, shole, transitions, param);
 
timing_reconstruction = toc(time_reconstruction);
 
timing_tot = toc(timetot);
 
timing = [timing, timing_transition, timing_reconstruction, timing_tot];
 
if nargout>3
    info.Gfull = Gfull;
    info.qvec = qvecfin;
    info.timing = timing;
    info.param = param;
end 
 
 
end
 
 
 


