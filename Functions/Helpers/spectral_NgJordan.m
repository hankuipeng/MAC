function [idx,flag,sumD] = spectral_NgJordan(X,K,varargin)

s = 0;
ds = 0;
nr = 0;
nn = 7;

if (rem(length(varargin),2) == 1)
	error('Optional parameters should always go by pairs');
else
	for i=1:2:(length(varargin)-1),
		if ~ischar (varargin{i}),
			error (['Unknown type of optional parameter name (parameter names must be strings).']);
		end
		
		switch lower(varargin{i})
		case 'sigma'
			s = double(varargin{i+1});
		case 'distance'
			ds = 1;
		case 'normalised'
			nr = 1;
		case 'nn'
			nn = double(varargin{i+1});
		otherwise
			error(['Unrecognized parameter: ''' varargin{i} '''']);
		end
	end
end

% Normalise if not distance matrix and not normalised
if ~nr && ~ds,
    % find columns that have no variation
    % and remove them
	active = max(X) - min(X);
	X = X(:, active~=0);
	
	% normalise in [-1,1]
	X = bsxfun(@minus,X, mean(X));
	X = X/max(max(abs(X)));
end

%% centralize and scale the data
%X = X - repmat(mean(X),size(X,1),1);
%X = X/max(max(abs(X)));

if ~ds,
	if s==0,
		% Zelnik Manor local scaling
		A = squareform(pdist(X));
		s = sort(A);
		% select distance to nn-th nearest neighbor
		s = s(nn+1,:);
		s(s<sqrt(eps)) = sqrt(eps);

		A = exp( -(A.^2)./(s'*s) );
		%assert( sum(diag(A) == zeros(n,1))==n);
	else
		A = squareform( exp(-(pdist(X)./(sqrt(2)*s)).^2) );
	end

% X is output of pdist function
elseif size(X,1)==1,
	if s==0,
		s = sort(squareform(X));
		% select distance to 7th nearest neighbor
		s = s(nn+1,:);
		s(s<sqrt(eps)) = sqrt(eps);

		A = exp( -(squareform(X).^2)./(s'*s) );
	else
		A = squareform( exp(-(X./(sqrt(2)*s)).^2) );
	end

% X is n \times n distance matrix
elseif size(X,1) == size(X,2) && issymmetric(X),
	if s==0,
		s = sort(X);
		% select distance to 7th nearest neighbor
		s = s(nn+1,:);
		s(s<sqrt(eps)) = sqrt(eps);

		A = exp( -(X.^2)./(s'*s) );
	else
		A = exp(-(X./(sqrt(2)*s)).^2);
	end
else 
	error('spectral_NgJordan: X matrix provided as input is incompatible\n');
end
% diagonal = 0
A(1:size(A,1)+1:end) = 0;

% to avoid division by zero cap degrees smaller than eps
deg = 1./sqrt( max([eps*ones(size(A,1),1), sum(A,2)], [], 2) );
%L = eye(n) -  (D*D').*A;

% select largest in magnitude eigenvalues
[U,lambdas,flag] = eigs((deg*deg').*A, K, 'lm');
clear A deg;

% Ng et al. 2002: Normalise rows to unit length
nrm = sqrt(sum(U.^2,2));
empty = find(nrm<sqrt(eps));
U(empty,:) = 1/sqrt(K);
nrm(empty) = 1.0;
U = bsxfun(@rdivide, U, nrm);

%fprintf('Eigen-Decomposition Complete!\n');
% k-means clustering
[idx,~,sumD] = kmeans(U,K,'Replicates',10);

