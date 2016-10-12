function [C,grpEvac,Cbin] = cum_evac(t_evac,grpIDs,groupProtocol,tbins,varargin)
% makes cumulative evacuation trajectories
% need to incorporate shleter space!!!!

%grpIDs = lowify(grpIDs,0);
[p,~] = size(grpIDs);
grpEvac = zeros(max(unique(grpIDs)),p).*NaN;
C = zeros(size(grpEvac)).*NaN;
for px=1:p  % for each sample / random grpID assignment
  for grp = removeval(unique(grpIDs)',0)  % for each group
    switch groupProtocol  % find evacuation # that causes group to go
        case 'ind',       tx_grp = find(grpIDs(px,:)==grp,1,'first');
        case 'individual',tx_grp = find(grpIDs(px,:)==grp,1,'first');
        case 'lastToGo',  tx_grp = find(grpIDs(px,:)==grp,1,'last');
        case 'lTG',       tx_grp = find(grpIDs(px,:)==grp,1,'last');
        case 'firstToGo', tx_grp = find(grpIDs(px,:)==grp,1,'first');
        case 'fTG',       tx_grp = find(grpIDs(px,:)==grp,1,'first');
        case 'majorityRule', 
            tx_grp = find(grpIDs(px,:)==grp);
            tx_grp = tx_grp(floor((numel(find(grpIDs(px,:)==grp))+1)/2));
        case 'mR',
            tx_grp = find(grpIDs(px,:)==grp);
            tx_grp = tx_grp(floor((numel(find(grpIDs(px,:)==grp))+1)/2));
    end
    if tx_grp <= size(t_evac,1)
        grpEvac(grp,px) = t_evac(tx_grp,px);
        if ~isnan(t_evac(tx_grp,p))
          C(grp,px) = sum(~isnan(grpEvac(:,px)));
        end
    end
  end
end
grpEvac = sort(grpEvac);
C = sort(C).*numel(find(grpIDs(px,:)==grp));

if nargin<2 %% nargout?
    Cbin = [];
else
    Cbin = histc(grpEvac,tbins);
    Cbin = cumsum(Cbin).*numel(find(grpIDs(px,:)==grp));
end

% add ss 
if any(ismemvar(varargin,'ss'))
    ss = varargin{find(ismemvar(varargin,'ss'),1,'first')+1};
    if ss<max(max(Cbin))
        Cbin(Cbin>ss) = ss;
    end
    if ss<max(max(C))
        C(C>ss) = NaN;
    end
end

end

