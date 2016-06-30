function [transitions, qvec] = ai_find_transitions(G,sthole,finhole,weight_disthole, weight_transition, weight_diffdist,nprop, exclude)
%AI_FIND_TRANSITIONS Find two optimal transitions
%   
%   Inputs parameters:
%       G           : Graph of transitions
%       sthole      : Starting node for the hole
%       finhole     : End node for the hole
%       weight_disthole : Regularization parameter - distance to the hole
%       weight transition : Regularization parameter - quality of the transitions
%       weight_diffdist : Regularization parameter - change in length
%       nprop       : Number of propositions
%       exclude     : Nodes to exclude 
%   Outputs parameters
%       transitions : Transitions
%       qvec        : Quality measurements of the transitions
%
%   This function searches for two ear-frendly transitions around the
%   hole. To do so, it solves an optimization problem.
%
 
% Authors: Nathanael Perraudin, Nicki Hollighaus
% Date   : June 2016
 
 
if nargin<4
    weight_disthole = 1;
end
 
if nargin<5
    weight_transition = 100;
end
 
if nargin<6
    weight_diffdist = 1;
end
 
if nargin<7
    nprop = 1;
end
 
if nargin<8
    exclude = [];
end
 
% Find starting and finishing positions
indices = find(G.d);
 
% All possible starting positions
ind_bfin = indices(indices<sthole);
% All possible finishing positions
ind_afin = indices(indices>finhole);
 
% We perform two jumps
[fin1,jump1,w1] = find(G.W(:,ind_bfin));
start1 = ind_bfin(jump1); % we need a reordering!
 
% find all possible jumps for the finishing nodes
[fin2,jump2,w2] = find(G.W(:,ind_afin));
start2 = ind_afin(jump2);
 
% Check if there is one empty set
if ~numel(w1) 
    disp('No path found! No connections before the hole.')
    transitions{1} = nan(2,2);
    qvec = [];
    return
end
if ~numel(w2) 
    disp('No path found! No connections after the hole.')
    transitions{1} = nan(2,2);
    qvec = [];
    return
end
 
 
 
% Let us now take care of the 3 penalizations
% 1) The distance to the hole
%    It is stored inside pmatrix.
p1 = sthole - ind_bfin;
p2 = ind_afin - finhole ;
%pmatrix = repmat(p1(start1),1,length(start2))' + repmat(p2(start2),1,length(start1));
pmatrix = repmat(p1(jump1),1,length(jump2))' + repmat(p2(jump2),1,length(jump1));
 
% 2) The weights of the edges
%    It is stored inside p2matrix.
 
% We would like a relative penalization, so we start by normalizing the
% weights
w1 = w1/max(w1);
w2 = w2/max(w2);
p2matrix = repmat(1./w1,1,length(w2))' + repmat(1./w2,1,length(w1));
 
% 3) The difference of the length of the two pieces
%    It is stored inside d
 
% The length of the first part
d1 = repmat(fin2,1,length(fin1)) - repmat(fin1,1,length(fin2))';
d1(d1<0) = nan; % Remove all negative lengths
% The length of the second part
d2 = repmat(start2,1,length(start1)) - repmat(start1,1,length(start2))';
if sum(sum(d2<0))
    error('Something is not normal!!!!')
end
d = abs(d1-d2); % Shall we square this term?
 
% remove fin1<sthole and fin2>finhole
% We do not want the hole to be in the middle of the recovery part
d((repmat(fin1,1,length(fin2))<=finhole)' & (repmat(fin2,1,length(fin1))>=sthole)) = nan;
 
% Remove the solutions that are overlapping with the excluded set
for ii = 1:size(exclude,1)
    d((repmat(fin1,1,length(fin2))<=exclude(ii,2) )' & (repmat(fin2,1,length(fin1))>=exclude(ii,1))) = nan;
end
 
% Objective function
obj = weight_diffdist*d + weight_disthole * pmatrix + weight_transition * p2matrix;
 
% Find the minimums
transitions = cell(nprop,1);
qvec = cell(nprop,1);
qvec{:} = nan(1,6);
for ii = 1:nprop
    [m1,i1] = min(obj);
    [testv,i2] = min(m1);
    
    qvec{ii} = [testv, ...
        weight_transition * p2matrix(i1(i2),i2), ...
        weight_diffdist * d(i1(i2),i2), ...
        weight_disthole * pmatrix(i1(i2),i2), ...
        weight_transition./w1(i2), ...
        weight_transition./w2(i1(i2))];
    if isnan(testv)
        if ii==1
            disp('Warning: The algorithm is not able to find a solution')
            transitions{ii} = nan(2,2);
            qvec = {nan(1,6)};
        else
            transitions = transitions(1:(ii-1));
            qvec = qvec(1:(ii-1));
            disp('Warning: The algorithm is not able to find so many solutions')           
        end
        break
    else
 
    transitions{ii}(1,1) = start1(i2);
    transitions{ii}(1,2) = fin1(i2);
    transitions{ii}(2,1) = fin2(i1(i2));
    transitions{ii}(2,2) = start2(i1(i2));
    
    % Remove the solution
    obj(i1,:) = nan;
    obj(:,i2) = nan;
    end
end
 
end


