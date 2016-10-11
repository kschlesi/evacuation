function par = index_each_cell(prange,pix)
% prange is an N-element cell array
% pix is a pxN matrix.
% par is a pxN matrix, where par(i,n) = prange{i}(pix(i))

if iscell(pix)
    ndim = length(pix);
    nindx = length(pix{1});
    pix = cell2mat(pix);
    pix = reshape(pix,[nindx,ndim]);
end
par = zeros(size(pix));
assert(numel(prange)==size(pix,2));
for p=1:size(pix,1)
  for i=1:numel(prange)
    par(p,i) = prange{i}(pix(p,i));
  end
end