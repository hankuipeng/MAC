%% load necessary paths
cd ..
repo_path = pwd;
addpath(genpath(repo_path))

clear repo_path


%% read in the data 
X = load('tfidf_full.csv');
Truth = load('num_labels.csv');


%% initial parameters
N = size(X, 1);
P = size(X, 2);
K = length(unique(Truth));


%% MAC
tic
grps_mac = MAC(X, K);
time_mac = toc;

cluster_performance(grps_mac, Truth)


%% save the results 
save('Experiments/MAC_results.mat')

