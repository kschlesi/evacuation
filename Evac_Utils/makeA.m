function [A,A_prefax] = makeA(N,qvals,ss)

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

% fix As for ss<N
if nargin>2
    if ss<N
      for i=1:qs
        % zero blocked transitions and transfer probabilities to allowed ones
        A(ss+1,:,i) = A(ss+1,:,i) + sum(A(ss+2:end,:,i));
        A(ss+2:end,:,i) = 0; 
        A(:,ss+2:end,i) = 0;
        % zero and re-compute diagonal
        A(:,:,i) = A(:,:,i).*(~eye(size(A(:,:,i)))); assert(all(~diag(A(:,:,i))));
        A(:,:,i) = A(:,:,i) - diag(sum(A(:,:,i))); %assert(all(~sum(A)));
      end
    end
end


end