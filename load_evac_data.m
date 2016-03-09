function [TotalPhitSeen,TotalPhitEvac,q,Crs_avg,missing,...
          Q1,T1,P1,C1,...
          Qm,Tm,Pm,Cm] = load_evac_data(makeFigs,indx,bins)
% This function loads evacuation data from 'evacuate_data.m'.
% Inputs:
%     makeFigs = boolean, whether to make figures
%     indx = optional indices of trials to include (default: individual, shelter=50)
%     bins = bin edges for Phit rounding (argument for 'histcounts', default [-0.05:0.1:1.05])
% Outputs:
%     TotalPhitSeen = total # of times Phits in each range were seen by
%                     subjects AT HOME. summed over all subjects and trials
%                     in 'indx'. (this is the H matrix)
%     TotalPhitEvac = total # of times any individual on any trial evacuated
%                     at Phits in each range. (this is the J matrix)
%     q = TotalPhitEvac./TotalPhitSeen. (this is Sean's 'theta' matrix)
%     Crs_avg = average cumulative # evacuated at each Phit value.
%     
% NOTES: loads newest versions of matrices ('evacuate_data.m') as of 12/14/15
%        uses fixed rounding bins (rounds to NEAREST 0.1)
%              (there are no Phits exactly equal to 0.x5)
%        makes sure evac times are appropriately indexed in gameinfo
%              (same Phit for t=1 and t=0, uses gameinfo matrix without t=0 value)
%        checks that no 'missing' trials are included
%        checks that evacTime and evacProb cross-reference correctly with gameinfo

% set time points, load data
tvec = 1:1:60;
load('evacuate_data.mat');
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
% import experiment PHIT Values
Q1 = gameinfo;
[ntrials, nts] = size(Q1); % row numbers = trial number, column numbers = time (1:1:60)
%%%%Q1(:,1) = 0.5;  % DO NOT manually fix initial Phit = 0.5 on all trials
assert(nts==length(tvec));

% import experiment Time at which individual evacuated
T1 = evacuateTime;
[N,~] = size(T1); % row numbers = individual number , column numbers = trial number
T1(T1==-1) = NaN; % replace -1 entries with NaN to indicate no evac
T1(:,all(T1==0)) = NaN; % replace 0's with NaN (in trials with no evacuation)

% import experiment Phit value at which the individual evacuated
P1 = evacuateProb;
P1(P1==-1) = NaN; % replace -1 entries with NaN to indicate no evac
P1(:,all(P1==0)) = NaN; % replace 0's with NaN (in trials with no evacuation)
try assert(sum(isnan(P1(:)))==sum(isnan(T1(:))));
catch  % fix instances where evac data NaNs in P1 and T1 do not match
    [i,j] = ind2sub(size(P1),find(isnan(P1(:))~=isnan(T1(:))));
    strlist = [];
    for ix=1:numel(i)
        strlist = sprintf([strlist,'\n\ttrial ',num2str(j(ix)),', subj ',num2str(i(ix)),...
                      ' at t = ',num2str(find(Q1(j(ix),:)==P1(i(ix),j(ix)),1,'first')),...
                      ', Phit = ',num2str(floor(P1(i(ix),j(ix))*10)/10)]);
    end
    warning([num2str(numel(i)),' inconsistent records: ',strlist,10,...
            'Treating as non-evacuations']);
    P1(i,j) = NaN;
    T1(i,j) = NaN;
    assert(sum(isnan(P1(:)))==sum(isnan(T1(:))));
end

% import experiment cumulative number evacauted
C1 = evapeocumu;

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% ****** WE WILL CONTANTINATE THE MATRICIES INTO TRIALS WE ******
% *************** ONLY WANT TO CONSIDER  *************************
% trials relevent to our scope: 16 trials total out of 160
% relevent trials are 19,28,29,49,64,67,76,102,108,112,113,123,125,130,143,144
if nargin<2
    % set default trial indices = individual, shelter=50 (shown above)
    indx = [19,28,29,49,64,67,76,102,108,112,113,123,125,130,143,144];
end
% check that all trials are valid and not in 'missing'
try assert(all(indx<=ntrials))
catch err
    error('Not enough trials to create requested matrix!');
end
try assert(~any(ismember(indx,missing)))
catch err
    warning('Missing trials requested!');
end

% ############################################################################################################%
% choosing only trials requested
Pm = P1(:,indx);
Tm = T1(:,indx);
Cm = C1(indx,:);
Qm = Q1(indx,:);

% we need to make 3 vectors, for each individual, we need to plot the
% number of times they saw a Phit value for the entire time they were in
% that trial, we need to plot the number of times they evacuated at a phit
% value in each trial, we need to plot the probability of each person to
% evacuate vs phit.

% This is how the PHIT Trajectroy was scaled during the experiment
if nargin<2
    %bins = 0:0.1:1.1; %% this is now incorrect! this rounds DOWN
    bins = -0.05:0.1:1.05; %% this rounds to NEAREST
end

% This is the total amount of times each individual seen a specified range
% of phit values.
% (commented-out sections are a check for evacs AFTER trial end; none found)
TotalPhitSeen = zeros(N,length(bins)-1);
% post_hit_evac = []; phe = false;
% phe_trial = []; phe_hit = [];
% phit hit is ranged from 0 to 1 in increments of 0.1
for i = 1:1:N % iterate over all individuals
    for j = 1:1:size(Qm,1) % iterate over relevant trials
        evacTime = Tm(i,j); is0 = 0;
        if isnan(evacTime) % no evacuation; ALL Phits seen until Qm = 1 or 0
            assert(Qm(j,end)==0||Qm(j,end)==1);
            h1 = histcounts(Qm(j,1:find(Qm(j,:)==Qm(j,end),1,'first')),bins);
%             if makeFigs
%               figure(1);hold on;plot(Qm(j,1:find(Qm(j,:)==Qm(j,end),1,'first')));hold off; title('no evac');
%             end
        else % evacuation; look only at Phits seen before evac decision
            if ~evacTime % if evacTime=0, use evacProb at t=1 (they are equal)
                evacTime = evacTime+1; is0 = 1;
            end
%             try assert(~(Qm(j,evacTime)==0||Qm(j,evacTime)==1));
%             catch  % instances where evac happened at Phit=1 or Phit=0
%                 warning(['evac at Phit = ' num2str(Qm(j,evacTime)) ...
%                          ', trial ' num2str(j) ', subj ' num2str(i)]);
%                 phe = true;
%                 post_hit_evac = [post_hit_evac evacTime-is0-find(Qm(j,:)==Qm(j,end),1,'first')];
%                 phe_trial = [phe_trial j];
%                 phe_hit = [phe_hit Qm(j,end)];
%             end
            assert(Qm(j,evacTime)==Pm(i,j));
            h1 = histcounts(Qm(j,1:evacTime),bins);
%             if makeFigs
%               figure(2);hold on;plot(Qm(j,1:evacTime+1));hold off; title('evac');
%             end
        end
        TotalPhitSeen(i,:) = TotalPhitSeen(i,:)+h1;
    end
end
% plot post-trial-end evacuations (none found)
% if phe && makeFigs
%     figure; hold on;
%     scatter(phe_trial(~phe_hit),post_hit_evac(~phe_hit),'ob');
%     scatter(phe_trial(~~phe_hit),post_hit_evac(~~phe_hit),'or');
%     title('# steps post-evac of phe decisions');
%     legend('no hit','hit','location','southoutside');
% end

% This is the total amount of times each individual evacuated in a
% specified range of phit values.
TotalPhitEvac = zeros(N,length(bins)-1);
for i = 1:1:N % iterate over all individials
    TotalPhitEvac(i,:) = TotalPhitEvac(i,:) + histcounts(Pm(i,:),bins);
end

% This is the probability for each individual to evacuate in the specified
% Phit range 
q = TotalPhitSeen.\TotalPhitEvac;

% ############################################################################################################
% compute the average cumulative number evacuated per Phit
[~,Qix] = histc(Qm,bins);
Crs_avg = zeros(size(bins));
for b=1:length(bins)
    Crs_avg(b) = mean(Cm(Qix==b));
end

% ############################################################################################################
% ############################################################################################################

% end of function modified by Kimberly





% we need to construct the total number of individual at home and evacaute
% for each PHIT

% Total number of participant who are ATHOME and observe PHIT.
% this should be equivalent to sum(TotalPhitSeen).
H = sum(TotalPhitSeen);

% Total number of participant who have Evacauted and observe PHIT.
% this should be equivalent to sum(TotalPhitEvac).
J = sum(TotalPhitEvac);

% Find time which experiment ended "Qm"
TrialEndTime = zeros(1,size(Qm,1)); % we will store the time at which the trials ended
for j = 1:size(Qm,1) % iterate through each trial
    TrialEndTime(j) = find(Qm(j,:)==Qm(j,end),1,'first');
end

% H = zeros(1,10); 
% J = zeros(1,10); 
% for i = 1:1:10 % iterate through the possible PHIT ranges {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1}
%     for j = 1:50 % iterate through each individual
%         for r = 1:16 % iterate through each trial
%             maxTime = TrialEndTime(r); % this is the time the trail (r) ended
%             IndvEvacTime = Tm(j,r); % this is the time the individual (j) evacuated in trial (r)
%             for k = 1:maxTime % iterate until trial (r) is done
%                 if k < IndvEvacTime && isnan(IndvEvacTime) ~= 1% iterate until the individual evacuates
%                     if i >= 1 && i < 10
%                         if Qm(r,k) < i/10 && Qm(r,k) >= (i-1)/10 % this makes sure we are in the correct PHIT range
%                             H(i) = H(i) + 1;
%                         end
%                     end
%                     if i == 10
%                         if Qm(r,k) <= i/10 && Qm(r,k) >= (i-1)/10 % this makes sure we are in the correct PHIT range
%                             H(i) = H(i) + 1;
%                         end
%                     end
%                        
%                 end
%                 if k >= IndvEvacTime && k <= maxTime % iterate until indv evacutes and trial (r) ends
%                     if i >= 1 && i < 10
%                         if Qm(r,k) < i/10 && Qm(r,k) >= (i-1)/10 % this makes sure we are in the correct PHIT range
%                             J(i) = J(i) + 1;
%                         end
%                     end
%                     if i == 10
%                         if Qm(r,k) <= i/10 && Qm(r,k) >= (i-1)/10 % this makes sure we are in the correct PHIT range
%                             J(i) = J(i) + 1;
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end

% ############################################################################################################
% ############################################################################################################

if makeFigs
    close all
%================================================================================================
% This is the CUMULATIVE NUMBER EVACAUTED vs PHIT for all
% trials relevent to our scope: 16 trials total out of 160
% relevent trials are 19,28,29,49,64,67,76,102,108,112,113,123,125,130,143,144

figure()
for i = [19,28,29,49,64,67,76,102,108,112,113,123,125,130,143,144]
    scatter(Q1(i,:), C1(i,:));
    hold on % this will hold each plot onto the same figure
    title(['Cumulative Number Evacauted for relevent trials ','  N = ', num2str(N-1)]);
    ylabel('Number Evacuated');
    xlabel('Phit');
    axis([0,1,0,50]);
end

% We will use this data to get the optimal parameters for the thoeretical
% model, then compare to the [AVERAGED CUMULATIVE NUMBER EVACUATED]
% parameters
%================================================================================================


%================================================================================================
% This is the AVERAGED CUMULATIVE NUMBER EVACUATED vs PHIT 

% We made use if the concatination, sorting and binning of Q1

Phit = 0:0.1:1; % Our binning had 11 bins, thus, increment by 0.1

figure()
scatter(Phit, Crs_avg)
title(['Averaged Cumulative Number Evacauted ', '  N = ', num2str(N)])
ylabel('Number Evacuated');
xlabel('Phit');
axis([0,1,0,50]);

% We will use this data to get the optimal parameters for the thoeretical
% model, then compare to the [CUMULATIVE NUMBER EVACUATED]
% parameters
%================================================================================================


%================================================================================================
% This is the PROBABILITY FOR EACH INDIVIDUAL TO EVACUATE vs PHIT, all in
% on plot

N = 50; % Total number of individuals
% create Phit range vector
Ph = 0:0.1:0.9;
% create Phit range matrix
PH = repmat(Ph,N,1);

figure()
for n = 1:1:N
    scatter(PH(n,:),q(n,:));
    hold on;
    title('\bf Probability for each Individual to Evacaute vs Phit');
    ylabel('Probability of Evacuation');
    xlabel('PHIT');
    axis([0,1,0,1]);
end

%================================================================================================


% %================================================================================================
% % This is the PROBABILITY FOR EACH INDIVIDUAL TO EVACUATE vs PHIT, All
% % individual plots
% 
% N = 50; % Total number of individuals
% a = floor(N/10); % # of Rows
% b = floor(N/(a*a)); % # of Columns
% Y = floor(N/a); % The # of individuals in each subplot figure
% % create Phit range vector
% Ph = 0:0.1:0.9;
% % create Phit range matrix
% PH = repmat(Ph,N,1);
% 
% % If we want to split plots into mutliple figures say each subplot figure
% % has Y plots, we need to make a for loop that has every Y plots in it
% 
% for j = 0:1:(a-1) % We want a total # of {a} figures
%     figure('OuterPosition',[200 0 1000 800]); % [left botton width height]
%     for n = 1:Y % Iterate over all Y # of individuals in each subplot
%         if n+j*Y > N % Stop plotting after the last individual is done
%             break;
%         else
%             h(n) = subplot(a,b,n); % make {a} rows / {b} columns 
%             subplot(h(n)); % n is the position for each subplot
%             scatter(PH(n+j*Y,:),q(n+j*Y,:));
% 
%             % create figure() legend for each subplot
%             leg = cell(1,1); % 1 - by - 1 cell for each subplot
%             for i=1:Y % Iterate over every subplot in the figure {we have Y # of plots in each figure}
%                 leg{i} = ['Indiv ' num2str(i + j*Y) ];
%             end
%             legend(leg(n),'Location','NorthWest');
%         end
%         % Make overall title for each subplot figure
%         ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%         text(0.5, 1,'\bf Probability of Evacuation vs PHIT, for Each Individual', 'FontSize', 20, 'HorizontalAlignment' ,'center','VerticalAlignment', 'top')
%     end
% end
% 
% %================================================================================================


%================================================================================================
% This is the CUMULATIVE NUMBER EVACAUTED vs TIME

figure()
plot(tvec, C1)
title(['Cumulative Number Evacauted ', 'Trial ', num2str(z), '  N = ', num2str(N-1)])
ylabel('Number Evacuated');
xlabel('Time (seconds)');
axis([tvec,1,0,50]);

% This data value will be compared to the mean number evacuated in the
% theoretical data
%================================================================================================


%================================================================================================
% This is PHIT vs TIME

figure();
plot(tvec, Q1);
title(['Phit Trajectory ', 'Trial ', num2str(z), '  N= ',num2str(N-1)]);
ylabel('Probability of Disaster to Strike');
xlabel('Time (seconds)');
axis([tvec,1,0,1]);

%================================================================================================
end % end makeFigs
 
% %================================================================================================
% % Gives rows of P on y-axis and columns of P on x-axis
% figure();
% bcolor(P)
% title(['Probability of Evacuation ', 'Trial ', num2str(z), '  N = ', num2str(N-1)]);
% ylabel('Time t');
% xlabel('Number Evacuated n');
% colorbar('location','eastoutside');
% 
% 
% 
% % Gives rows of P on y-axis but inverted and columns of P on x-axis
% figure();
% imagesc(0:N, T_range, P)
% title(['Probability of Evacuation ', 'Trial ', num2str(z), '  N = ', num2str(N-1)]);
% ylabel('Time t');
% xlabel('Number Evacuated n');
% colorbar('location','eastoutside');
% %================================================================================================


end  % end main function
