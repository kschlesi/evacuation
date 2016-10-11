function [sortmat,ix] = sortentry(inmat,dir,which,sortvec)
% given a 2D matrix 'inmat,' this function treats either rows or columns as
% linked entries (indicated by 'dir') and sorts the entries by the value of
% a given field 'which.' Default: columns are entries; 'which' names a row
% by whose values to sort. Returns: sorted matrix 'sortmat' (same size as
% 'inmat') and vector 'ix' giving the rearranged index order of the entries
% Optional 4th entry 'sortvec': set 'which' = 0

if strcmp(dir,'row')
    inmat = inmat';
end

% go to sortvec mode?
if nargin>3 % need numel(sortvec) to be size(inmat,1)
    if ~which && (sum((size(sortvec)>1))==1) && (numel(sortvec)==size(inmat,1))
        [useless,ix] = sort(sortvec);
        if size(ix,1)==1
            ix = ix';
        end
        which = 1;
    else
        error('Set sorting field = 0 or reshape sorting vector!');
    end
else
% if NO sortvec
    [useless,ix] = sort(inmat,1);
end

% create 'sortmat' from ix (sorted 'which' field or input vector)
sortmat = zeros(size(inmat));
for i=1:size(inmat,1)
    sortmat(i,:) = inmat(ix(i,which),:);
end

ix = ix(:,which);

if strcmp(dir,'row')
    sortmat = sortmat';
end

end