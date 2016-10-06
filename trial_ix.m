function [ix,ss,gs] = trial_ix(groupProtocol,shelterSpace,groupSize,missing)

load('experiment.mat');

% special case: if one argument structured [gp,ss,gs] then return str for
% protocol and other arguments as looked up in gamestatus info
if nargin==1
    gp = gamestatus(groupProtocol,1);
    ss = gamestatus(groupProtocol,2);
    gs = gamestatus(groupProtocol,3);
    switch gp
        case 0, ix = 'individual';
        case 1, ix = 'firstToGo';
        case 2, ix = 'lastToGo';
        case 3, ix = 'majorityRule';
    end
    return
end

% normal thing: 
switch groupProtocol
    case 'ind',
        if nargin<3
            groupSize = 1;
        end
        gp = 0;
    case 'individual', 
        if nargin<3
            groupSize = 1;
        end
        gp = 0;
    case 'firstToGo', gp = 1;
    case 'fTG', gp = 1;
    case 'lastToGo', gp = 2;
    case 'lTG', gp = 2;
    case 'majorityRule', gp = 2;
    case 'mR', gp = 3;
end
ix1 = find(gamestatus(:,1)==gp);
ix2 = find(gamestatus(:,2)==shelterSpace);
ix3 = find(gamestatus(:,3)==groupSize);
ix = intersect(intersect(ix1,ix2),ix3);
ss = []; gs = [];

if nargin==4
    ix = removeval(ix,missing);
end

end