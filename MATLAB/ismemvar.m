function [Lia,Locb] = ismemvar(A,B,varargin)
% this function works like 'ismember,' but allows operation on cell arrays 
% with both string and non-string cells. 
% *****
% NOTE: if a combination of cell array + non-cell array is passed in that 
% would cause an error in 'ismember' (i.e., the cell array is not a cell 
% array of strings OR the non-cell array is not a string), the DEFAULT
% behavior of this function is to check whether the contents of each cell 
% in the cell array are identical to the ENTIRE non-cell array input.
% EXAMPLES: ismemvar({'one'; 'two'; 1; 2}, 'one') --> [1; 0; 0; 0]
%           ismemvar({'one'; 'two'; [1 2]}, [1 2]) --> [0; 0; 1]
%           ismemvar({'one'; 'two'; 1; 2}, [1 2]) --> [0; 0; 0; 0]

% cases in which it directly mimics 'ismember'
if (~iscell(A) && ~iscell(B)) || ...
    ((ischar(A) || iscellstr(A)) && (ischar(B) || iscellstr(B)))
    [Lia,Locb] = ismember(A,B,varargin{:});
else
    % otherwise.... use cellmem. Locb is now supported.
    if ~iscell(A); A = {A}; end;
    if ~iscell(B); B = {B}; end;
    [Lia,Locb] = cellmem(A,B);
    %Locb = [];
    if numel(varargin)
      warning('rows or legacy behavior not supported for these arguments'); 
    end
end

end