function bns = bar_unpaired(ncat,varargin)
% this function plots a combined bar histogram of various pools of values
% that are not necessarily of the same size. first arg, 'ncat,' gives
% number of pools of values to histogram (each pool in a different color).
% first ncat entries in 'varargin' should be vectors or matrices of values
% to hist. other optional argument pairs: 'nbins', scalar number of bins.

% load passed-in value of nbins or set default=10
if any(ismemvar(varargin,'nbins'))
    nbins = varargin{find(ismemvar(varargin,'nbins'),1,'first')+1};
else
    nbins = 10;
end

% load all values to hist
reshape_row = @(in_) reshape(in_,1,numel(in_));
allVals = cell2mat(cellfun(reshape_row,varargin(1:ncat),'UniformOutput',false))';

% find hist bins and sort into those bins
[~,bns] = hist(allVals,nbins);
hist_bins = @(in_) hist(in_(:),bns)';
plotMat = cell2mat(cellfun(hist_bins,varargin(1:ncat),'UniformOutput',false));

% plot bin counts as bars
bar(plotMat);
set(gca,'XTickLabel',bns);

end