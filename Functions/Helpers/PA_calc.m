function [theta, diff, dmax] = PA_calc(F, G, prop)

% Note (2018/11/03):
% This function calculates subspace similarity based on the proposed
% subspace similarity measure in MAC. The principal angles are calculated
% using the method proposed in 'Principal angles between 
% subspaces in an a-based scalar product'.

% prop: the proportion of variability in the data to be expalained
% maximally by the PCs.

% Last updated: 1 Oct 2019


%% F-related quantities
if size(F,1) == 1
    Qf=F';
    df = size(Qf,2);
else
    [U1 S1 V1] = svd(F','econ'); % U1 contains the eigenvectors
    vals1 = diag(S1).^2;
    k1 = sum(cumsum(vals1)/sum(vals1)<=prop); 
    if k1==0
        k1=1;
    end
    Qf = U1(:,1:k1);
    df = size(Qf,2);
end



%% G-related quantities
if size(G,1)==1
    Qg=G';
    dg = size(Qg,2);
else
    [U2 S2 V2] = svd(G','econ'); % U1 contains the eigenvectors
    vals2 = diag(S2*S2');
    k2 = sum(cumsum(vals2)/sum(vals2)<=prop); 
    if k2==0
        k2=1;
    end
    Qg = U2(:,1:k2);
    dg = size(Qg,2);
end

if dg > df
    % swap matrices
    cp = Qf;
    Qf = Qg;
    Qg = cp;
    
    % swap dims
    cp = df;
    df = dg;
    dg = cp;
end

%% reduced SVD  
[Y S Z] = svd(Qf'*Qg, 'econ');

%% final step
theta = acos(diag(S));
diff = abs(df-dg);
dmax = max(df,dg);
u = Qf*Y;
v = Qg*Z;

end