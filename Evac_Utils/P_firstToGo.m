function prob = P_group(n,ng,N,Gs)
% calculates the probability, in a system with N individuals in groups of
% size Gs, that n individuals selected randomly without replacement will be
% chosen from exactly ng different groups. this probability is P(n,ng).

if any(mod([n,ng,N,Gs],1))
    error('all inputs must be integers!');
end
if mod(N,Gs)
    error('does not support groups that do not divide N evenly!');
end

% ng must be a viable number of groups
if ng > N/Gs
    prob = 0;
    return
end

% base case: P(n,1)
if ng == 1
    i = 1:1:n-1;
    prob = prod((Gs-n+i)./(N-n+i));
    return
end

% base case: P(n,n)
if n == ng
    i = 0:
end




end