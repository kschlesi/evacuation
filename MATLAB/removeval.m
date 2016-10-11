function shortlist = removeval(A,val)
% a function that removes all entries of the specified value from vector A
% returns 'shortlist', a vector of length (length(A) - times 'val' appears in A)

if ~ismatrix(A) || ~( numel(A)==size(A,1) || numel(A)==size(A,2) )
    error('Input must be a row or column vector!')
end

if ~ismatrix(val) || ~( numel(val)==size(val,1) || numel(val)==size(val,2) )
    error('Input values in a row or column vector!')
end

shortlist = A;
for i=1:numel(val)
    shortlist = shortlist(shortlist~=val(i));
    %disp(shortlist);
end

end