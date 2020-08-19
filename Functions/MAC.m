%%% Inputs:
% X: N by P data matrix.
% K: number of clusters.
% prop: the proportion of variability to capture in a subspace with the
% basis vectors.
% thres: the tolerance parameter in RREF
% normalise: 1 or 0, depending on whether to normalise the data or not.

% Last updated: 19th Aug. 2020


function finlabels = MAC(X, K, prop, thres, normalise)

if nargin < 3 || isempty(prop)
    prop = 0.8;
end

if nargin < 4
    thres = 0;
end

if nargin < 5
    normalise = 1;
end


%% normalise the data to unit length 
if normalise == 1
    X0 = X;
    X = norml2(X0, 1);    
end


%% get the dimensionalities
N = size(X, 1);
P = size(X, 2);


%% Obtain the Reduced Row Echelon Form (RREF)
Xrref = rref(X');
Xrref = decimal(Xrref, 4); % keep the precision up to 4 decimal points


%% assign those data objects that are connected through RREF together 
Adj = zeros(N,N);

for i = 1:(N-1)
    for j = (i+1):N
        if Xrref(:,i)'*Xrref(:,j) ~= thres
            Adj(i,j) = 1;
        end
    end
end

Adjacency = Adj + Adj';
 
comp = conncomp(graph(Adjacency), 'OutputForm', 'cell');


%% create one sub-matrix for each current cluster
submat = {};

for kk = 1:size(comp,2)
    submat{kk} = X(comp{kk},:);
end


%% Part 3: create a matrix of dissimilarities
Dfull = PA_dis(submat, P, prop);
groups = SC_Jordan(Dfull, K);


%% Part 4: recover the full labels
finlabels = zeros(N, 1);

for k = 1:K
    gind = find(groups == k);
    sz = size(gind,1);
    for s = 1:sz
        finlabels(comp{gind(s)}) = k;
    end
end


end
