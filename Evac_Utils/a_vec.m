function [avec,aix] = a_vec(t_,P_,A_,Phit_traj_)

amat = A_(:,:,floor(Phit_traj_(floor(t_)+1)*10+1)).*repmat(P_',length(P_),1);
aix = tril(ones(length(P_)),-1);
avec = amat(~~aix);

end