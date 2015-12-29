function Ptrans = grp_trans(N,Gs,groupProtocol)
% this returns an (Ng+1) x (N+1) matrix, where entry (i,j) gives 
% P_group(j,i,N,Gs) for i=0:Ng (Ng=N/Gs) and j=0:N

Ng = N/Gs;
Ptrans = zeros(Ng+1,N+1);

if strcmp(groupProtocol,'fTG') || strcmp(groupProtocol,'firstToGo')
    for j=0:Ng
        for i=0:N
            jx=j+1;
            ix=i+1;
            Ptrans(jx,ix) = P_group(i,j,N,Gs);
        end
    end
else
    Ptrans = [];
end

end