function ix = find_each_cell(carray,vals,num,pos)
% crray is an N-element cell array
% vals is an N-element vector.
% ix is an (num x N) matrix, where ix(i,n) = find(carray{i}==vals(i),num,pos)
    % default num = 1, default pos = 'first'
    
if nargin<4
    pos = 'first';
end
if nargin<3
    num = 1;
end
N = numel(carray);
ix = zeros(num,N);
for i=1:N
    ix(:,i) = find(carray{i}==vals(i),num,pos);
end


end