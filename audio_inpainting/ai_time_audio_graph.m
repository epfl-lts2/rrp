function [G, Gfull,param, timing] = ai_time_audio_graph(s,fs,param)
%AI_TIME_AUDIO_GRAPH Create a graph of similarities from an audio signal
%   Usage: G = ai_time_audio_graph(s, fs)
%          G = ai_time_audio_graph(s, fs, param)
%          [G, Gfull] = ai_time_audio_graph(...)
%          [G, Gfull, param] = ai_time_audio_graph(...)
%          [G, Gfull, param, timing] = ai_time_audio_graph(...)
%
%   Input parameters:
%       s           : signal
%       fs          : sampling frequency
%       param       : structure of optional parameters
%   Output parameters:
%       G           : graph of audio transitions
%       Gfull       : a graph with more audio transitions
%       param       : structure of optional parameters (updated)
%       timing      : timing (features computation, graph construction)
%
%   This function create a graph of similarity inside a song with the
%   following steps:
%
%   1. Downsampling of the song 
%   2. Compute audio features from a Time-Frequency transform
%   3. Create the graph of nearest neighbors
%   4. Select only relevant connections
%
%   *param* is a structure of optional parameter containing the following
%   fields:
%
%   * *param.fsmax* : maximum sampling frequency (default 12 000 Hz). The
%     algorithm will cut everything above this frequency before any other
%     operation.
%
%   * *param.win_length* : window length (default: 1024)
%   * *param.M* : number of frequency channels (default: *param.win_length*)
%   * *param.a* : shift in time (default: 128)
%   * *param.win* : window used (default: 'itersine')
%
%   * *param.dbrange* : range of dB used in the spectrograph (default: 50)
%
%   * *param.use_flann* : use flann library for the computation of the
%     graph (default 1)
%   * *param.k*         : number of neighbors for the graph (default 10)
%   * *param.loc*       : locality parameter (default 0). How much being
%     far away from each other is important?
%
%   * *param.hsize*     : numbeer of time bin for one patch (default 4).
%
%   * *param.diagdist*     : size of the convolution kernel  (default 20).
%   * *param.threshold*    : threshold to cut small values (default 2).
%   * *param.featuretype*  : type of features (default 0), see the function
%     ai_compute_features for more details.
%
%   References: perraudin2016audio
%
 
% Authors: Nathanael Perraudin, Nicki Hollighaus
% Date   : June 2016
 
 
% 0) Handling optional parameters
if nargin < 3
    param = struct;
end
 
if ~isfield(param,'fsmax'), param.fsmax = 12000; end
if ~isfield(param,'verbose'), param.verbose = 0; end
 
 
if ~isfield(param,'use_flann'), param.use_flann = 1 ; end
if ~isfield(param,'k'), param.k = 40 ; end
if ~isfield(param,'graph_type'), param.graph_type = 'knn' ; end
if ~isfield(param,'loc'), param.loc = 0 ; end
 
if ~isfield(param,'hsize'), param.hsize = 0 ; end
 
if ~isfield(param,'diagdist'), param.diagdist = 40 ; end
if ~isfield(param,'threshold'), param.threshold = 2 ; end
if ~isfield(param,'featuretype'), param.featuretype = 4 ; end
 
if isfield(param,'subgraph')
    
    if ~isfield(param.subgraph,'time_in')
        fprintf('Missing starting time. Use the full graph \n');
        use_subgraph = 0;
    elseif ~isfield(param.subgraph,'time_fin')
        fprintf('Missing final time. Use the full graph \n');
        use_subgraph = o;
    else
        if ~isfield(param.subgraph,'time_range')
            param.subgraph.time_range = 5;
        end
    use_subgraph = 1;
    end
    
else
    use_subgraph = 0;
end
 
 
%% STEP 1
% Downsampling , we do not really care about high frequencies + making
% the sound mono (we can use stereo sound later on)
if size(s,2)>1
    s = (s(:,1)+s(:,2))/2;
end
 
% % Original size of the signal
% Noriginal = length(s);
 
% Downsample if necessary 
if fs > param.fsmax
    ratio = ceil(fs/param.fsmax);
    s = resample(s,1,ratio);
    fs = fs/ratio;
end
 
 
%% STEP 2
% 2) Create audio features
if param.verbose
    fprintf('Compute local audio features\n');
end
% featuremat = sim_features(s,fs,param);
% 
% %% Create the feature patches
% 
% 1) Compute useful quantities
 
 
% if nargout > 3
%     [patches,patchmax] = sim_features(s,fs,param);
% else 
%     patches = sim_features(s,fs,param);
% end
tfeatures = tic;
 
[featuremat, param, patchmax] = ai_compute_features(s,fs,param);
 
hsize = param.hsize;
marg = floor(hsize / 2);  % left - right margin
 
if hsize
    [vsize,owidth] = size(featuremat);
    dim = hsize*vsize;
 
    % 2) Initialize variables    
    patches = zeros(dim+1,owidth-2*marg);
 
    for jj = 1 : owidth-2*marg    
        patches(1:end-1,jj) = reshape(...
            featuremat(:, (1:hsize)+(jj-1) ),[],1);
        % for  distance of the coordinates
        patches(end,jj) = jj*param.loc;
    end
else
    patches = featuremat;
end
timing_features = toc(tfeatures);
 
%% Step 3
% Create the time variable
N = size(patches,2);
% Check this !!
time = ((marg):(marg+N-1) ) * param.a * ratio;
 
%% STEP 4
tgraph = tic;
% Create the graph of nearest neighbors
gparam.use_flann = param.use_flann;
gparam.type = 'knn';
gparam.k = max(param.k,floor(N/1000)); % Patch graph minimum number of connections (KNN).
 
gparam.center = 0;
gparam.rescale = 0;
gparam.resize = 0;
gparam.symetrize_type = 'full';
gparam.type = param.graph_type;
gparam.target_degree = param.k;
 
if use_subgraph
    t_ini = max(round(param.subgraph.time_in/ratio/param.a - param.subgraph.time_range*fs/param.a  ),1);
    t_inf = round(param.subgraph.time_in/ratio/param.a );
    t_outf = min(round(param.subgraph.time_fin/ratio/param.a + param.subgraph.time_range*fs/ param.a),N);
    t_outi = round(param.subgraph.time_fin/ratio/param.a);
    if param.verbose
        fprintf('Compute graph with %d vertices\n' , size(patches,2));
    end
    Gd = gsp_nn_part_graph(patches',[1:t_inf,t_outi:size(patches,2)],[t_ini:t_inf,t_outi:t_outf],gparam);
else
    if param.verbose
        fprintf('Compute graph with %d vertices\n' , size(patches,2));
    end
    Gd = gsp_nn_graph(patches', gparam);    
end
    % Update plotting parameters
    x = sin(2*(1:Gd.N)/Gd.N*pi*0.95);    
    y = cos(2*(1:Gd.N)/Gd.N*pi*0.95);
    Gd.coords = [x;y]';
    Gd.plotting.limits = [-1.1 1.1 -1.1 1.1];
    Gd = gsp_graph_default_parameters(Gd);
 
%% STEP 5
%  Select only relevant connections. To do so we create another graph G2 !
 
if param.verbose
    fprintf('Select only good connections\n');
end
% 1) remove the diagonal
W = Gd.W;
% W( logical( 1-(triu(ones(Gd.N),param.diagdist) + tril(ones(Gd.N),-param.diagdist))) ) = 0;
% W( logical(spdiags(ones(Gd.N,param.diagdist*2+1), -param.diagdist:param.diagdist, Gd.N, Gd.N)) ) = 0;
W( logical(spdiags(ones(Gd.N,2*round(param.diagdist/2)+1), ...
    -round(param.diagdist/2):round(param.diagdist/2), Gd.N, Gd.N)) )...
    = 0;
 
    
    
% 2) We want to have multiple connections in time to assert a good
% edge. So we compute a convolution with a diagonal kernel. This will
% diffuse the energy of isolate connections and keep strong multiple
% connections. 
dd = param.diagdist;
convk = diag(tripuls(-round(dd)/2:round(dd)/2,round(dd)+1)); % kernel
 
if use_subgraph
     Wconv = sconv2(W,convk,'same');
%     Wconv = zeros(size(W));
%     Wconv(:,[t_ini:t_inf,t_outi:t_outf]) = conv2(full(W(:,[t_ini:t_inf,t_outi:t_outf])),convk,'same');
else
    %Wconv = conv2(full(W),convk,'same');
    Wconv = sconv2(W,convk,'same');
 
end
 
% Remove small elements
[row,col,v] = find(Wconv);
row = row(v>param.threshold);
col = col(v>param.threshold);
v = v(v>param.threshold);
Wconv = sparse(row,col,v,Gd.N,Gd.N);
 
% figure(1); imagesc(log(Wconv))
 
% Remove the starting and final connections. They are wrong because the
% signal is not periodic
Nr = round(param.win_length/param.a);
Wconv(1:Nr,:) = 0;
Wconv(:,1:Nr) = 0;
Wconv((end-Nr):end,:) = 0;
Wconv(:,(end-Nr):end) = 0;
 
 
if nargout>1
    % Add new parameter to the graph structure
    Gfull = Gd;
    Gfull.W = Wconv;
    Gfull = gsp_graph_default_parameters(Gfull);
    Gfull.time = time;
    Gfull.ratio = ratio;
    Gfull.patchmax = patchmax;
    Gfull.fs = fs;
    Gfull.a = param.a;
    Gfull.features = featuremat;
 
end
 
 
 
 
if param.verbose>1
    figure(102)
    imagesc(W)
    colorbar
    title('Nearest Neighboor weight matrix')
    figure(103)
    imagesc(Wconv)
    colorbar
   % colormap(flipud(gray))
    title('Weight matrix after convolution')
end
 
 
% 3) Compute local maxima for this
C1 = Wconv(2:(end-1),2:(end-1));
C2 = Wconv(1:(end-2),3:(end));
C3 = Wconv(3:(end),1:(end-2));
C4 = ((C1-C2)>0) & ((C1-C3)>0);
 
Wconv(2:(end-1),2:(end-1)) = Wconv(2:(end-1),2:(end-1)).* C4;
Wconv(1,:) = 0;
Wconv(:,1) = 0;
Wconv(end,:) = 0;
Wconv(:,end) = 0;
 
 
 
if param.verbose>1
    figure(104)
    imagesc(Wconv)
    colorbar
    title('Nearest Neighboor weight after first maxima selection')
end
 
 
C1 = Wconv(2:(end-1),2:(end-1));
C2 = Wconv(1:(end-2),1:(end-2));
C3 = Wconv(3:(end),3:(end));
C4 = ((C1-C2)>0) & ((C1-C3)>0);
 
Wconv(2:(end-1),2:(end-1)) = Wconv(2:(end-1),2:(end-1)).* C4;
Wconv(1,:) = 0;
Wconv(:,1) = 0;
Wconv(end,:) = 0;
Wconv(:,end) = 0;
 
if param.verbose>1
    figure(105)
    imagesc(Wconv)
    colorbar
    title('Nearest Neighboor weight after second maxima selection')
end
% 4) Create the new graph
G = Gd;
G.W = Wconv;
G = gsp_graph_default_parameters(G);
clear Gd
 
% Add new parameter to the graph structure
G.time = time;
G.ratio = ratio;
G.patchmax = patchmax;
G.fs = fs;
G.a = param.a;
G.features = featuremat;
 
timing_graph = toc(tgraph);
 
timing = [timing_features, timing_graph];
end

