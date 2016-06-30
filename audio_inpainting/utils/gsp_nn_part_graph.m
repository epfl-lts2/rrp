function [ G ] = gsp_nn_part_graph( X,ind1,ind2, param )
%GSP_NN_PART_GRAPH Create a nearest neighbors subgraph from a point cloud
%   Usage :  G = gsp_nn_part_graph( X,ind1,ind2 );
%            G = gsp_nn_part_graph( X,ind1,ind2, param );
%
%   Input parameters:
%       X           : Input points
%       ind1        : Indices 1
%       ind2        : Indices 2
%       param       : Structure of optional parameters
%
%   Output parameters:
%       G           : Resulting graph
%
%   'gsp_nn_part_graph( X,ind1,ind2, param )' creates a subgraph from
%   positional data. The points are connected to their neighbors (either
%   belonging to the k nearest neighbors or to the epsilon-closest
%   neighbors. This function only looks for connections between the nodes
%   indiced by *ind1* and *ind2*.
%
%   Additional parameters
%   ---------------------
%
%   * *param.type*      : ['knn', 'radius']   the type of graph (default 'knn')
%   * *param.use_flann* : [0, 1]              use the FLANN library
%   * *param.center*    : [0, 1]              center the data
%   * *param.rescale*   : [0, 1]              rescale the data (in a 1-ball)
%   * *param.sigma*     : float               the variance of the distance kernel
%   * *param.k*         : int                 number of neighbors for knn
%   * *param.epsilon*   : float               the radius for the range search
%   * *param.use_l1*    : [0, 1]              use the l1 distance
%
%   See also: gsp_nn_graph gsp_nn_distanz
%

% Author: Johan Paratte, Nathanael Perraudin
% Date: 16 June 2014
% Testing: test_rmse

    if nargin < 4
    % Define parameters
        param = {};
    end
    
    %Parameters
    if ~isfield(param, 'type'), param.type = 'knn'; end
    if ~isfield(param, 'use_flann'), param.use_flann = 0; end
    if ~isfield(param, 'center'), param.center = 1; end
    if ~isfield(param, 'rescale'), param.rescale = 1; end
    if ~isfield(param, 'k'), param.k = 10; end
    if ~isfield(param, 'epsilon'), param.epsilon = 0.01; end
    if ~isfield(param, 'use_l1'), param.use_l1 = 0; end
    if ~isfield(param, 'target_degree'), param.target_degree = 0; end;

  
    paramnn = param;
    paramnn.k = param.k +1;
    [indx, indy, dist, Xout, ~, epsilon] = gsp_nn_distanz(X(ind1,:)',X(ind2,:)',paramnn);
    Xout = transpose(Xout);
    
    
    switch param.type
        case 'knn'
            if param.use_l1
                if ~isfield(param, 'sigma'), param.sigma = mean(dist); end
            else
                if ~isfield(param, 'sigma'), param.sigma = mean(dist)^2; end
            end
        case 'radius'
            if param.use_l1
                if ~isfield(param, 'sigma'), param.sigma = epsilon/2; end
            else
                if ~isfield(param, 'sigma'), param.sigma = epsilon.^2/2; end
            end
        otherwise
            error('Unknown graph type')
    end
    
    n = size(X,1);
    
    if param.use_l1
        W = sparse(ind1(indx), ind2(indy), double(exp(-dist/param.sigma)), n, n);
    else
        W = sparse(ind1(indx), ind2(indy), double(exp(-dist.^2/param.sigma)), n, n);
    end
    
    % We need zero diagonal
    W(1:(n+1):end) = 0;     % W = W-diag(diag(W));
  
    
    
    
    % Computes the average degree when using the epsilon-based neighborhood
    if (strcmp(param.type,'radius'))
        text = sprintf('Average number of connections = %d', nnz(W)/size(W, 1));
        disp(text);
    end

    % Precopy the necessary values
    
	W = gsp_symmetrize(W,'full');   
    %Fill in the graph structure
    G.N = n;
    G.W = W;
    G.coords = X;
    %G.limits=[-1e-4,1.01*max(x),-1e-4,1.01*max(y)];
    if param.use_l1
        G.type = 'Part nearest neighbors l1';
    else
        G.type = 'Part nearest neighbors';
    end
    %G.vertex_size=30;
    G.sigma = param.sigma;
    G = gsp_graph_default_parameters(G);
end









% 
% 
% 
% 
% 
% 
% 
% function [ G ] = gsp_nn_part_graph( X,ind1,ind2, param )
% %GSP_NN_PART_GRAPH Create a nearest neighbors graph from a point cloud
% %   Usage :  G = gsp_nn_part_graph( X,ind1,ind2 );
% %            G = gsp_nn_part_graph( X,ind1,ind2, param );
% %
% %   Input parameters:
% %       X           : Input points
% %       ind1        : Indices 1
% %       ind2        : Indices 2
% %       param       : Structure of optional parameters
% %
% %   Output parameters:
% %       G           : Resulting graph
% %
% %   'gsp_nn_part_graph( X,ind1,ind2, param )' creates a subgraph from
% %   positional data. The points are connected to their neighbors (either
% %   belonging to the k nearest neighbors or to the epsilon-closest
% %   neighbors. This function only looks for connections between the nodes
% %   indiced by *ind1* and *ind2*.
% %
% %   Additional parameters
% %   ---------------------
% %
% %   * *param.type*      : ['knn', 'radius']   the type of graph (default 'knn')
% %   * *param.use_flann* : [0, 1]              use the FLANN library
% %   * *param.center*    : [0, 1]              center the data
% %   * *param.rescale*   : [0, 1]              rescale the data (in a 1-ball)
% %   * *param.sigma*     : float               the variance of the distance kernel
% %   * *param.k*         : int                 number of neighbors for knn
% %   * *param.epsilon*   : float               the radius for the range search
% %   * *param.use_l1*    : [0, 1]              use the l1 distance
% %
% %   See also: gsp_pointcloud
% %
% 
% % Author: Johan Paratte, Nathanael Perraudin
% % Date: 16 June 2014
% % Testing: test_rmse
% 
%     if nargin < 2
%     % Define parameters
%         param = {};
%     end
%     
%     %Parameters
%     if ~isfield(param, 'type'), param.type = 'knn'; end
%     if ~isfield(param, 'use_flann'), param.use_flann = 0; end
%     if ~isfield(param, 'center'), param.center = 1; end
%     if ~isfield(param, 'rescale'), param.rescale = 1; end
%     if ~isfield(param, 'k'), param.k = 10; end
%     if ~isfield(param, 'epsilon'), param.epsilon = 0.01; end
%     if ~isfield(param, 'use_l1'), param.use_l1 = 0; end
%     if ~isfield(param, 'target_degree'), param.target_degree = 0; end;
% 
%         
%     
%     % test if the binaries of flann are working
%     if param.use_flann
%        try
%             paramsflann.algorithm = 'kdtree';
%             paramsflann.checks = 32;
%             paramsflann.trees = 1;
%             tmp = rand(100,10);
%             [NN, D] = flann_search(tmp', tmp', 3, paramsflann);
%        catch
%             warning('Flann not compiled, going for the slow algorithm!')
%             param.use_flann = 0;
%        end
%     end
%     
%     Xin = X;
%     [Nel, Nfeatures] = size(Xin);
%     Nel2 = numel(ind2);
%     Xout = Xin;
%     
%     if strcmpi(param.type,'knn') && param.k> Nel-2
%         error('Not enougth samples')
%     end
%     
%     %Center the point cloud
%     if (param.center)
%         Xout = Xin - repmat(mean(Xin), [Nel, 1]);
%     end
%     
%     
%     
%     %Rescale the point cloud
%     if (param.rescale)
%         bounding_radius = 0.5 * norm(max(Xout) - min(Xout));
%         scale = nthroot(Nel, min(Nfeatures, 3)) / 10;
%         Xout = Xout .* (scale / bounding_radius);
%     end
%     
%     switch param.type
%         %Connect the k NN
%         case 'knn'
%             k = param.k;
%      
%             spi = zeros(Nel2*k,1);
%             spj = zeros(Nel2*k,1);
%             spv = zeros(Nel2*k,1);
%             
%             %Find kNN for each point in X (Using a kdtree)
%             if param.use_flann
%                 if param.use_l1
%                     error('Not implemented yet')
%                 end
%                 paramsflann.algorithm = 'kdtree';
%                 %TODO : optimize parameters in function of the number of
%                 %points
%                 paramsflann.checks = 256;
%                 paramsflann.trees = 1;
%                 % Use flann library
%                 [NN, D] = flann_search(Xout(ind1,:)', Xout(ind2,:)', k, paramsflann);
%                 NN = transpose(NN);
%                 D = transpose(D);
%             else
%                 %Built in matlab knn search
%                 if ~isreal(Xout)                   
%                     Xout2 = [real(Xout),imag(Xout)];
%                     if param.use_l1
%                         kdt = KDTreeSearcher(Xout2(ind1,:), 'distance', 'cityblock');
%                         [NN, D] = knnsearch(kdt, Xout2(ind2,:), 'k', k , 'Distance','cityblock');
%                     else
%                         kdt = KDTreeSearcher(Xout2(ind1,:), 'distance', 'euclidean');
%                         [NN, D] = knnsearch(kdt, Xout2(ind2,:), 'k', k );
%                     end                   
%                 else
%                     if param.use_l1
%                         kdt = KDTreeSearcher(Xout(ind1,:), 'distance', 'cityblock');
%                         [NN, D] = knnsearch(kdt, Xout(ind2,:), 'k', k ,  'Distance','cityblock');
%                     else
%                         kdt = KDTreeSearcher(Xout(ind1,:), 'distance', 'euclidean');
%                         [NN, D] = knnsearch(kdt, Xout(ind2,:), 'k', k );
%                     end                   
%                 end
%             end
%             
% 
%             if param.use_l1
%                 if ~isfield(param, 'sigma'), param.sigma = mean(D(:)); end
%             else
%                 if ~isfield(param, 'sigma'), param.sigma = mean(D(:))^2; end
%             end
%             
%             % Fill the 3-col values with [i, j, exp(-d(i,j)^2 / sigma)]
%             for ii = 1:Nel2
%                 spi((ii-1)*k+1:ii*k) = repmat(ind2(ii), k, 1);
%                 spj((ii-1)*k+1:ii*k) = ind1(NN(ii, :))';
%                 if param.use_l1
%                     spv((ii-1)*k+1:ii*k) = exp(-D(ii,:)/param.sigma);
%                 else
%                     spv((ii-1)*k+1:ii*k) = exp(-D(ii,:).^2/param.sigma);
%                 end
%             end
%             
% 
%             
%         %Connect all the epsilon-closest NN
%         case 'radius'
%             %Create KDTree for fast NN computation
%             
%             if param.use_l1
%                 kdt = KDTreeSearcher(Xout(ind1,:), 'distance', 'cityblock');
%             else
%                 kdt = KDTreeSearcher(Xout(ind1,:), 'distance', 'euclidean');
%             end
%             
%             if (param.target_degree == 0) 
%                 epsilon = param.epsilon;
%             else
%                 target_d = floor(param.target_degree);
%                 [NN, D] = knnsearch(kdt, Xout(ind2,:), 'k', target_d);
%                 avg_d2 = mean(D(:,end).^2);
%                 epsilon = sqrt(avg_d2);
%             end
%             tic;
%             %Find all neighbors at distance <= epsilon for each point in X
%             if param.use_l1
%                [NN, D] = rangesearch(kdt, Xout(ind2,:), epsilon,'Distance' ,'cityblock');
%             else
%                [NN, D] = rangesearch(kdt, Xout(ind2,:), epsilon, 'distance', 'euclidean' );
%             end
%             toc;
%             
%             %Counting non-zero elements
%             count = 0;
%             for ii = 1:Nel2
%                count = count + length(NN{ii}) - 1; 
%             end
%             spi = zeros(count,1);
%             spj = zeros(count,1);
%             spv = zeros(count,1);
%             start = 1;
%             
% 
%             % Fill the 3-col values with [i, j, exp(-d(i,j)^2 / sigma)]
%             for ii = 1:Nel2
%                 len = length(NN{ii});
%                 spi(start:start+len-1) = repmat(ind2(ii), len, 1);
%                 spj(start:start+len-1) = ind1(NN{ii}');
%                 if param.use_l1
%                     if ~isfield(param, 'sigma'), param.sigma = epsilon/2; end
%                     spv(start:start+len-1) = exp(-D{ii}/param.sigma);
%                 else
%                     if ~isfield(param, 'sigma'), param.sigma = epsilon^2/2; end
%                     spv(start:start+len-1) = exp(-D{ii}.^2/param.sigma);
%                 end
%                 start = start + len;
%             end
%             
%         otherwise
%             error('Unknown type : allowed values are knn, radius');
%     end
% 
%     %Actually create the sparse matrix from the 3-col values
%     W = sparse(spi, spj, spv, Nel, Nel);
%     
%     % Computes the average degree when using the epsilon-based neighborhood
%     if (strcmp(param.type,'radius'))
%         text = sprintf('Average number of connections = %d', nnz(W)/size(W, 1));
%         disp(text);
%     end
% 
%     % Precopy the necessary values
%     
% 	W = gsp_symetrize(W,'full');   
%     %Fill in the graph structure
%     G.N = Nel;
%     G.W = W;
%     G.coords = Xout;
%     %G.limits=[-1e-4,1.01*max(x),-1e-4,1.01*max(y)];
%     if param.use_l1
%         G.type = 'Part nearest neighbors l1';
%     else
%         G.type = 'Part nearest neighbors';
%     end
%     %G.vertex_size=30;
%     G.sigma = param.sigma;
%     G = gsp_graph_default_parameters(G);
% end
% 