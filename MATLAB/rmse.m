function err = rmse(tru,pred)
% given a set of true measures and a set of predictions of those measures,
% this function computes the root mean squared error between the two sets.
% if more than one non-singleton dimension, this takes the mean along the
% first dimension (e.g., averaging over each column if two dimensions).
% if only one imput, the function computes rmse with predictions of 0.

if nargin==1
    pred = zeros(size(tru));
end

try assert(numel(tru)==numel(pred))
catch
    error('RMSE inputs must have the same number of elements');
end

if length(tru)==numel(tru)
    tru = tru(:);
end

if length(pred)==numel(pred)
    pred = pred(:);
end

err = sqrt(nanmean(((tru-pred).^2)));

end