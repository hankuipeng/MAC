function finlabels = MAC_stage2(X, Xrref, K, thres, prop)

N = size(X, 1);
P = size(X, 2);

if nargin < 5 || isempty(prop)
    prop = 0.8;
end


%% assign those data objects that are connected through RREF together 
Adj = zeros(N,N);

for i = 1:(N-1)
    for j = (i+1):N
        if Xrref(:,i)'*Xrref(:,j) > thres
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