function Q = newmansQ(atest,Ctest)
% computes Newman's Q-value for modularity of a network with nxn adjacency
% matrix A, given a partition C.

P = full(sparse(Ctest,1:length(Ctest),1));
Q = sum(sum((P*atest).*P))/sum(sum(atest));

end