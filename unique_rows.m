function [umat,ix] = unique_rows(mat,sameNan)
% this function eliminates all non-unique rows in the 2D matrix 'mat'
% it returns ix, the indices of the (first appearances of) unique rows, and
%            umat, which is equal to mat(ix,:)
% if second argument 'sameNan' = 1, it treats all NaN as equivalent

if nargin<2
    sameNan = 0;
end

if sameNan
    nanpad = max(unique(mat))+1;
    mat(isnan(mat)) = nanpad;
end

ix = [];
for i=1:size(mat,1)-1
    ix = [ix; find(all(bsxfun(@eq,mat(i,:),mat(i+1:end,:)),2)) + i];
end
ix = removeval(1:1:size(mat,1),unique(ix));
umat = mat(ix,:);

if sameNan
    umat(umat==nanpad) = NaN;
end