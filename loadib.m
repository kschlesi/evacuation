function ib = loadib(subjects,missings)
% subjects = scalar number of subjects scanned
% missings = vector of subjects whose scans to exclude

% set defaults for choking and officer
    if nargin<2 || ~numel(missings)
        ib = [];
        if subjects == 22  % choking default
            ib = [1;3;4;6;7;8;9;10;11;12;13;14;15;16;17;19;20;21;22]; 
        end
        if subjects == 86  % officer default (hybrid fMRI)
            load('officer_subject_names.mat');
            assert(numel(ib)==86);
        end        
    else
        ib = removeval((1:subjects)',missings);
    end

end