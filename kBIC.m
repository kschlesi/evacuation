function BIC = kBIC(X, clust, C)
% inputs : X, an NxM data set (points x factors)
%          clust, an Nx1 kmeans clustering solution
%          C, a KxM matrix of means (centroids) of each cluster
% output : BIC, bayesian information criterion for model
%               assuming spherical Gaussian clusters

[N,M] = size(X); % number of points, number of dimensions
K = numel(unique(clust)); % number of clusters
Nk = zeros(K,1);
for k=1:K
    Nk(k) = sum(clust==k);
end
assert(sum(Nk)==N);
sigma = sum(sum((X - C(clust,:)).^2))./(M*(N-K));
BIC = sum(Nk.*log(Nk)) - N*log(N) - N*M*log(2*pi*sigma^2)/2 - M*(N-K)/2;

end