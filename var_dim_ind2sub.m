function pix = ind2sub_var_dim(sizemat,indx)

clear pix;
[pix{1:numel(prange)}] = ind2sub(sizemat,indx);
pix = cell2mat(pix);