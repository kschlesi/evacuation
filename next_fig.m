function nextNo = next_fig()
% determine next figure number not already used in an open figure

% get all open figs
hAllFigs = findobj(0,'Type','figure');

% sort figs numberically
if numel(hAllFigs)>1
    hAllNums = cell2mat(get(hAllFigs,'Number')); % Sort the figure Numbers
else
    hAllNums = get(hAllFigs,'Number');
end

% set next number
if numel(hAllNums);
    nextNo = min(removeval(1:max(hAllNums),hAllNums));
    if ~numel(nextNo)
        nextNo = max(hAllNums)+1;
    end
else
    nextNo = 1; 
end


end
