function prob = P_group(n,ng,N,Gs)
% calculates the probability, in a system with N individuals in groups of
% size Gs, that n individuals selected randomly without replacement will be
% chosen from exactly ng different groups. this probability is P(n,ng).
% NOTE: this works but it scales as 2^N with recursion as-is. FIX

if any(mod([n,ng,N,Gs],1))
    error('all inputs must be integers!');
end
if mod(N,Gs)
    error('does not support groups that do not divide N evenly!');
end

% ng must be a viable number of groups
if ng > N/Gs || ng<=0
    prob = 0;
    return
end

% limits on n
if n<=0 || n<ng || n>N
    prob = 0;
    return
end

% base case: P(n,1)
if ng == 1
    if n==1
        prob = 1;
        return
    end
    i = 1:1:n-1;
    prob = prod((Gs-n+i)./(N-n+i));
    return
end

% base case: P(n,n)
if n == ng
    i = 2:n;
    prob = prod((N-Gs*(i-1))./(N-i+1));
    return
end

% recursion case: P(n,n-i)
i = n-ng;
disp(['Calculating P(' num2str(n) ',' num2str(i) ')']);
prob = ( P_group(n-1,n-i,N,Gs)*((n-i)*Gs-n+1) + ...
         P_group(n-1,n-i-1,N,Gs)*(N-(n-i-1)*Gs) )/(N-n+1);

end