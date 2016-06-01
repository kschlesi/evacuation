function [emb_matrix,effdim] = CSE(A)
%CALC_CSEMBEDDING performs constant shift embedding
%                 of a similarity matrix A, which contains similarity
%                 measures between 0 and 1 for all samples
% 
% the INPUT is
%   A: nxn similarity matrix 
% the OUTPUT is
%   emb_matrix: (n-1)xn embedding matrix in (n-1) nominal dimensions
%               the column vectors are the vector values of each sample
%   effdim:     effective dimension (sum of eigenvalues)
%               in practice you might want to consider only the first
%               effdim components of the vector for each sample,
%               effectively taking emb_matrix(1:effdim,:)

% transforming similarities into dissimilarities
%diss_matrix=A-1;
diss_matrix = 1-A-eye(size(A));
try assert(all(diag(diss_matrix)==0));
catch err
    warning('Forcing all self-dissimilarities to 0');
    diss_matrix(~~eye(size(A))) = 0;
    assert(all(diag(diss_matrix)==0));
end
% determining sample size
[n,~]=size(diss_matrix);
% Constant Shift Embedding into an Euclidean Space of dimension n-1
Q=eye(n)-ones(n)/n;
Dt=diss_matrix-2*min(real(eig(-Q*diss_matrix*Q/2)))*(ones(n)-eye(n));
Sct=-Q*Dt*Q/2;
% calculation of the transformation eigenvalues and eigenvectors
[vec,val]=eig(Sct);
vec=real(vec);
val=real(val);
[dval,srt]=sort(diag(val),'descend');

p=n-1;
%p=n-1; THIS IS THE UPPER BOUND
%p=sum(dval>1e-4);
vec=vec(:,srt);
valp=diag(dval(1:p));
vecp=vec(:,1:p);
% calculation of the embedding matrix
emb_matrix=(vecp*(valp.^(1/2)))';
emb_matrix=real(emb_matrix);
% estimating the "effective dimension" (sum of eigenvalues)
effdim=real(trace(valp)/valp(1,1));
% now each column vector of emb_matrix determines the vector values of the
% corresponding sample from A
end