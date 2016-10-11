function [t_evac,grpIDs] = group_evac_times(N,p,groupSize,t_mat,r_mat,ix,jx)
% given N (# participants), p (# simulations), groupSize, and r_mat (matrix
% of simulated reaction IDs), this function returns grpIDs (a pxN matrix of
% randomized group IDs for the N participants) and t_evac (matrix of sorted
% times of evacuation for each 

    grpID1 = repmat(1:floor(N/groupSize),1,groupSize);
    grpID1 = [grpID1 zeros(1,N-length(grpID1))];
    grpIDs = zeros(p,N);
    n_evac = zeros(size(r_mat)); 
    n_evac(~~r_mat) = ix(r_mat(~~r_mat))-jx(r_mat(~~r_mat));
    t_evac = zeros(max(sum(n_evac)),p).*NaN;
    for px=1:p
      tx = [];
      for ne = removeval(unique(n_evac),0)'
         tx = [tx; repmat(t_mat(n_evac(:,px)==ne,px),ne,1)];
      end
      t_evac(1:sum(n_evac(:,px)),px) = sort(tx)';
      grpIDs(px,:) = grpID1(randperm(N));
    end