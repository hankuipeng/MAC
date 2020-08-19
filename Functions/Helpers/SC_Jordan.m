function groups = SC_Jordan(Dfull, K)

if size(Dfull,1) == K
    groups = 1:K;
else
    %[groups,~] = SpectralClustering(W, K);
    [groups,~,~] = spectral_NgJordan(Dfull, K);
    %[groups,~,~] = spectral_NgJordan(Dfull, K,'distance',1);
end

end