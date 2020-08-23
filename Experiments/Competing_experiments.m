% You will need the corresponding code to run the following algorithms,
% most of which are publicly available.


%% load MAC path
cd ..
repo_path = pwd;
addpath(genpath(repo_path))

clear repo_path


%% read in the data 
X = load('tfidf_full.csv');
Truth = load('num_labels.csv');
X_unit = norml2(X, 1);


%% initial parameters
N = size(X, 1);
P = size(X, 2);
K = length(unique(Truth));


%% LDA 
grps_lda = load('Experiments/grps_lda.csv');
cluster_performance(grps_lda, Truth)


%% KSC
% step 1: determine q
% step 2: KSCq


%% retrieve partial tfidf for each category
tfidf1 = X(Truth == 1, :);
tfidf2 = X(Truth == 2, :);
tfidf3 = X(Truth == 3, :);
tfidf4 = X(Truth == 4, :);
tfidf5 = X(Truth == 5, :);


%% eigen-decomposition of each cluster 
[Vec1, val1] = eig(cov(tfidf1)); Val1 = diag(val1); vals1 = sort(Val1, 'descend');
[Vec2, val2] = eig(cov(tfidf2)); Val2 = diag(val2); vals2 = sort(Val2, 'descend');
[Vec3, val3] = eig(cov(tfidf3)); Val3 = diag(val3); vals3 = sort(Val3, 'descend');
[Vec4, val4] = eig(cov(tfidf4)); Val4 = diag(val4); vals4 = sort(Val4, 'descend');
[Vec5, val5] = eig(cov(tfidf5)); Val5 = diag(val5); vals5 = sort(Val5, 'descend');


%% plot the eigenvalues in descending order
plot(vals1(1:500))
hold on 
plot(vals2(1:500))
hold on 
plot(vals3(1:500))
hold on 
plot(vals4(1:500))
hold on 
plot(vals5(1:500))
hold on 
legend({'Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Cluster 5'}, 'Location', 'northeast')
hold off


%% KSC
ksc_path = '/home/hankui/Dropbox/Clustering/WorkLog/By_Algorithm/SCAL/SCAL_git';
addpath(genpath(ksc_path))
clear ksc_path

max_iter = 100;
q = 200;

tic;
[grps_ksc, tot_re_q] = KSCq(X, K, max_iter, q);
time_ksc = toc;
cluster_performance(grps_ksc, Truth)


%% SSC 
ssc_path = '/home/hankui/Dropbox/Clustering/WorkLog/By_Algorithm/SSC_ADMM';
addpath(genpath(ssc_path));
clear ssc_path

tic;
r = 0; affine = 0; alpha = 20; outlier = 1; rho = 1;
[missrate, grps_ssc, CKSym] = SSC(X', r, affine, alpha, outlier, rho, Truth);
time_ssc = toc;
cluster_performance(grps_ssc, Truth)


%% LRR
lrr_path = '/home/hankui/Dropbox/Clustering/WorkLog/By_Algorithm/LRR';
addpath(genpath(lrr_path));
clear lrr_path

% grps_lrr = lrr_all(X', K, lambda);
% time_lrr = toc;
% cluster_performance(grps_lrr, Truth)

% obtain the coefficient matrix Z
tic;
lambda = 0.01;
Z = solve_lrr(X', X', lambda);
time_lrr1 = toc;

tic;
A = abs(Z) + abs(Z');

% spectral clustering 
grps_lrr = SpectralClustering(A, K);
time_lrr2 = toc;

cluster_performance(grps_lrr, Truth)
time_lrr = time_lrr1 + time_lrr2;


%% PDDP 
pddp_path = '/home/hankui/Dropbox/Clustering/WorkLog/By_Algorithm/PDDP';
addpath(genpath(pddp_path));
clear pddp_path

tic;
[grps_pddp, tree] = pddp(X, K);
time_pddp = toc;

cluster_performance(grps_pddp, Truth)

%plot(tree, Truth, lines(K))


%% SC(X)
tic;
[grps_scx] = spectral_NgJordan(X, K);
time_scx = toc;

cluster_performance(grps_scx, Truth)


%% SC(A)
tic;

% normalise the data to unit length 
X_unit = norml2(X, 1);    


% get the dimensionalities
N = size(X_unit, 1);
P = size(X_unit, 2);

% Obtain the Reduced Row Echelon Form (RREF)
[Q, R] = qr(X');
Xrref0 = rref(R);
Xrref = denoise_c(Xrref0);

% assign those data objects that are connected through RREF together 
Adj = zeros(N,N);

for i = 1:(N-1)
    for j = (i+1):N
        if Xrref(:,i)'*Xrref(:,j) ~= 0
            Adj(i,j) = 1;
        end
    end
end

A = Adj + Adj';

% spectral clustering
[grps_sca] = spectral_NgJordan(A, K, 'distance', 0);

time_sca = toc;

cluster_performance(grps_sca, Truth)


%%
