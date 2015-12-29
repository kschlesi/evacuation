function [A,A_prefax] = makeA(N,qvals)

qvals = qvals(:)';
qs = numel(qvals);
A_prefax = zeros(N+1,N+1);
A = zeros(N+1,N+1,qs);
for nx=1:N+1
    for ix=1:N+1
        n = nx-1;
        i = ix-1;
        if i < n  % constrain only to evacuation moves
            A_prefax(nx,ix) = nchoosek(N-i,n-i);
            A(nx,ix,:) = shiftdim((qvals.^(n-i)).*((1-qvals).^(N-n)),-1);
        end
    end
end
for j=1:qs
    A(:,:,j) = A(:,:,j).*A_prefax; assert(all(~diag(A(:,:,j))));
    A(:,:,j) = A(:,:,j) - diag(sum(A(:,:,j)));
end

end