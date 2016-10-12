function str = trial_conv(trialNo,ss,groupSize,groupProtocol)
% converts trial info into string following standard naming convention
% ignores (the fact that there are any) Missings right now.

if strcmp(groupProtocol,'Ind') || strcmp(groupProtocol,'Individual')
    ig = 'Ind';
    groupProtocol = [];
    groupSize = [];
    news = [];
else
    ig = 'Grp';
    news = '-';
end

trialNoNew = (find(trial_ix(groupProtocol,ss,groupSize)==trialNo));

if ~ischar(ss)
    ss = num2str(ss);
end

switch groupProtocol
    case 'fTG', gps = 'FTG';
    case 'firstToGo', gps = 'FTG';
    case 'lastToGo', gps = 'LTG';
    case 'lTG', gps = 'LTG';
    case 'mR',  gps = 'MV';
    case 'majorityRule',  gps = 'MV';
    case 'Individual', gps = 'Ind';
    
end

str = ['Trial ' ig ss '-' gps num2str(groupSize) news num2str(trialNoNew)];

end