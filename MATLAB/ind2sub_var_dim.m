function pix = ind2sub_var_dim(sizemat,indx,isCell)
% calls ind2sub, catches ndim outputs where ndim = length(sizemat), returns
% a numel(indx) x ndim matrix of these outputs

clear pix;
ndim = length(sizemat);
[pix{1:ndim}] = ind2sub(sizemat,indx);
if ~isCell || nargin<3
  pix = cell2mat(pix);
  pix = reshape(pix,[numel(indx),ndim]);
end

end