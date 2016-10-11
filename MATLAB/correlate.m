function [rr,pp] = correlate(sample1,sample2,varargin)

% convert samples to column vectors and check corresponding lengths
sample1 = sample1(:);
sample2 = sample2(:);
if (length(sample1)~=length(sample2))
    error('samples must have the same number of elements');
end

% set type of correlation coefficient (default = Pearson's r)
% other types: 'Spearman' or 'Kendall'
if any(ismemvar(varargin,'type'))
    type = varargin{find(ismemvar(varargin,'type'),1,'first')+1};
else
    type = 'Pearson';
end
% set symbol for chosen type
switch type
    case 'Pearson',  symbol = 'r';
    case 'Spearman', symbol = '\rho';
    case 'Kendall',  symbol = '\tau';
    otherwise,       symbol = [];
end

% compute correlation & p-value
if any(ismemvar(varargin,'partial'))
    partialout = varargin{find(ismemvar(varargin,'partial'),1,'first')+1};
    partialout = partialout(:);
    assert(length(partialout)==length(sample1));
    [rr,pp] = partialcorr(sample1,sample2,partialout,'type',type);
%     part_tag = varargin{find(ismemvar(varargin,'partial'),1,'first')+2};
else
    [rr,pp] = corr(sample1,sample2,'type',type);
%     part_tag = [];
end
% use the following lines for corrcoef only
% rr = rr(1,2);
% pp = pp(1,2);

% create scatterplot with correlation
if ~any(ismemvar(varargin,'NoFigure'))
    if ~any(ismemvar(varargin,'ExtFigureCmd'))
    figure;
    end
    if any(ismemvar(varargin,'DispResiduals'))
       resType = varargin{find(ismemvar(varargin,'DispResiduals'),1,'first')+1};
       lin1 = fitlm(partialout,sample1); 
       lin2 = fitlm(partialout,sample2);
       eval(['sample1 = lin1.Residuals.' resType ';']);
       eval(['sample2 = lin2.Residuals.' resType ';']);    
    end
    scatter(sample1,sample2);
    if ~any(ismemvar(varargin,'RemoveLeg'))
    legend([type '''s ' symbol ' = ' num2str(rr) ', p = ' num2str(pp)]);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    % otherwise.... use cellmem. Locb is NOT supported here.
    if ~iscell(A); A = {A}; end;
    if ~iscell(B); B = {B}; end;
    Lia = cellmem(A,B);
    Locb = [];
    if numel(varargin)
      warning('rows or legacy behavior not supported for these arguments'); 
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function comp = cellmem(A,B)
% this takes two cell arrays, A and B. returns a logical vector C the same
% size as A, with C(i) = 1 if the contents of A(i) are identical to the
% contents of ANY cell in B.

comp = zeros(size(A));
for a=1:numel(A)
    xinb = zeros(size(B));
    for b=1:numel(B)
        % if either cell (A or B) CONTAINS a char, use strcmp;
        % else use all(equals)
        if ischar(A{a}) || ischar(B{b})
            xinb(b) = strcmp(A{a},B{b});
        else
            try AeqB = A{a}==B{b};
            catch oops
                if oops; AeqB = 0; end;
            end
            xinb(b) = (all(AeqB(:)) && numel(A{a})==numel(B{b}));
        end
    end
    comp(a) = any(xinb(:));
end
        
end