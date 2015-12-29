function PbXStochastic_sample_data()
% Considering Experimental Data

% OBJECTIVE: We will plot cumulative vs phit, then do it the same way
%            but now concatinate it, then bin it, from that, we will look at
%            the averaged values. With this, we will costruct the linear
%            parameter (Slope), Threshold parameter (Mean Phit, Max
%            Evacuated), Hill parameters (Mean Phit, Max Evacuted,
%            Steepness of the S-Curve) for both the averaged and the
%            non-averaged plot

% MODIFICATIONS: We will look for optimal parameters now by limiting the
%                range of the best fit. Use maximum likelihood estimation
%                for this. Then with the thoery let q = linear, threshold,
%                Hill-function with respective parameters.

format long

% Time Variable: define range of timesteps in which to give solutions
T_range1 = [0, 60];   % TIme evolution is seconds of trial
T_range = 0:1:T_range1(2);

% MATLAB does not start at zero, so since we have 50 individual subject we
% need to make our final stop to be 51 
N = 51;  % number of total states 

%=========================================================================
% GENERAL MASTER EQUATION
% dP[n,t+dt] = SUM[i~=n:N](P[i,t]*P[n,t+dt|i,t] - P[n,t]*P[i,t+dt|n,t])

% If n > i then P[i, t + dt|n, t] = 0 "we cannot leave the shelter"
% If i > n then P[n, t + dt|i, t] = 0 "we cannot leave the shelter"

% we need P[n|i] = [N-i  n-i] q^(n-i) (1-q)^(N-n)
% where [N-i n-i] = (N-i)! / ((n-i)! * (N-n)!)
%=========================================================================

% ******Phit is a random stochastic process******
% make probability of one person to evacuate still linear to Phit but Phit
% is now a stochastic process that is given from the experiment data, so we
% need to import the file probability.csv

z = 125; % this is the trial number

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
% import experiment PHIT Values
Q1 = importdata('gameinfo.xlsx'); 
length(Q1);
size(Q1); % row numbers = trial number, column numbers = time (0:1:60)
Q = Q1(z,:); % this is a row vector 
length(Q);  % length(Q1) = total number of time-steps (0 to 60)
size(Q);

% we need to replace (-1) in evacuateTime by the time which the disaster
% hit or did not hit. Let i be the time in seconds at which the outcome of
% the event ended. We start off at 2 since at 1, value of Phit = Q(1) = 0
for i = 2:61
    if Q(i) == 1;
        Q(i); 
        i; % Time when disaster struck
        break;
    end
    if Q(i) == 0;
        Q(i); 
        i; % Time when disaster missed
        break;
    end
end
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% import experiment Time at which individual evacuated
T1 = importdata('evacuateTime.xlsx');
length(T1);
size(T1); % row numbers = individual number , column numbers = trial number
T = T1(:,z); % this is a column vector
length(T);  % length(T1) = total number of states "N" (0 to 50)
size(T);

% replace each element of T with value (-1) by NaN seconds
% !!!!!!!!! we then only make the axis show 0 to 60 seconds that way we can
% neglect this term, there might be a better way to approach this

% this will save the index of vector T at which the change will be made
X = []; 
for j = 1:length(T)
    if T(j) == -1;
        X = [X j]; % this will let me know which individual did not evacuate
        T(j) = NaN; % replace T(j) = -1 with NaN since Individual did not evacuate
    end
end
T; % this is the correct Time of evacuation vector

% this will sort the data of T in ascending order
[Ts,I] = sort(T);
length(Ts);
Ts; % this is the sorted T values
I;  % this is the corresponding indicies of T,
    % we will use these indicies to sort P since the indicies in the
    % unsorted form of P corresponds to the values in the unsorted form of
    % T
    
% % This is the sorted values of T with respect to the indicies of P
% Ts = T(I);
% length(Ts);
% Ts;
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    
% import experiment Phit value at which the individual evacuated
P1 = importdata('evacuateProb.xlsx');
length(P1);
size(P1); % row numbers = individual number , column numbers = trial number
P = P1(:,z); % this is a column vector
length(P);  % length(P1) = total number of states "N" (0 to 50)
size(P);

% replace each element of P with value (-1) by the value 1.1
% !!!!!!!!! we then only make the axis show 0 to 1, that way we can
% neglect this term, there might be a better way to approach this

% this will save the index of vector P at which the change will be made
Y = []; 
for j = 1:length(P)
    if P(j) == -1;
        Y = [Y j]; % this will let me know which individual did not evacuate
        P(j) = NaN; % replace P(j) = -1 with NaN since individual did not evacaute
    end
end
P; % this is the correct Probability of evacuation vector

% % this will sort the data of P in ascending order
% [Ps,I] = sort(P);
% length(Ps);
% Ps; % this is the sorted P values
% I;  % this is the corresponding indicies of P,
%     % we will use these indicies to sort T since the indicies in the
%     % unsorted form of T corresponds to the values in the unsorted form of
%     % P

% This is the sorted values of P with respect to the indicies of T
Ps = P(I);
length(Ps);
Ps;
% % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% import experiment cumulative number evacauted
C1 = importdata('evapeocumu.xlsx');
length(C1);
size(C1); % row numbers = individual number , column numbers = trial number
C = C1(z,:); % this is a row vector
length(C); % length(C1) = total number of time-steps (0 to 60)
size(C);
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


% ****** WE WILL CONTANTINATE THE MATRICIES INTO TRIALS WE ******
% *************** ONLY WANT TO CONSIDER  *************************
% trials relevent to our scope: 16 trials total out of 160
% relevent trials are 19,28,29,49,64,67,76,102,108,112,113,123,125,130,143,144

% ############################################################################################################
% ############################################################################################################
Y = zeros(50,160); 
for i = 1:1:50 % iterate through each individual
    for j = 1:1:160 % iterate through each trial
        if P1(i,j) == -1;
            Y(i,j) = i; % this will let me know which individual (i) did not evacuate in trial (j)
            P1(i,j) = NaN; % replace P(i,j) = -1 with NaN since individual did not evacaute
        end
    end
end
Y; % this is the index of matrix P1 for the individuals who did not evacuate
size(Y); % [50 - by - 160]
length(Y);

X = zeros(50,160); 
for i = 1:1:50 % iterate through each individual
    for j = 1:1:160 % iterate through each trial
        if T1(i,j) == -1;
            X(i,j) = i; % this will let me know which individual (i) did not evacuate in trial (j)
            T1(i,j) = NaN; % replace T(i,j) = -1 with NaN since individual did not evacaute
        end
    end
end
X; % this is the index of matrix T1 for the individuals who did not evacuate
size(X); % [50 - by - 160]
length(X);

% Concatinate matrix where we only pick specific columns and get matrix out
Pm = [P1(:,19) P1(:,28) P1(:,29) P1(:,49) P1(:,64) P1(:,67) P1(:,76) P1(:,102) P1(:,108)...
    P1(:,112) P1(:,113) P1(:,123) P1(:,125) P1(:,130) P1(:,143) P1(:,144)];
Tm = [T1(:,19) T1(:,28) T1(:,29) T1(:,49) T1(:,64) T1(:,67) T1(:,76) T1(:,102) T1(:,108)...
    T1(:,112) T1(:,113) T1(:,123) T1(:,125) T1(:,130) T1(:,143) T1(:,144)];
% Concatinate matrix where we only pick specific row and get matrix out
Cm = [C1(19,:);C1(28,:);C1(29,:);C1(49,:);C1(64,:);C1(67,:);C1(76,:);C1(102,:);C1(108,:);...
    C1(112,:);C1(113,:);C1(123,:);C1(125,:);C1(130,:);C1(143,:);C1(144,:)]; 
Qm = [Q1(19,:);Q1(28,:);Q1(29,:);Q1(49,:);Q1(64,:);Q1(67,:);Q1(76,:);Q1(102,:);Q1(108,:);...
    Q1(112,:);Q1(113,:);Q1(123,:);Q1(125,:);Q1(130,:);Q1(143,:);Q1(144,:)];

% we need to make 3 vectors, for each individual, we need to plot the
% number of times they saw a Phit value for the entire time they were in
% that trial, we need to plot the number of times they evacuated at a phit
% value in each trial, we need to plot the probability of each person to
% evacuate vs phit.

% This is how the PHIT Trajectroy was scaled during the experiment
% box 1:  0 <= PHIT < 0.1
% box 2:  0.1 <= PHIT < 0.2
% box 3:  0.2 <= PHIT < 0.3
% box 4:  0.3 <= PHIT < 0.4
% box 5:  0.4 <= PHIT < 0.5
% box 6:  0.5 <= PHIT < 0.6
% box 7:  0.6 <= PHIT < 0.7
% box 8:  0.7 <= PHIT < 0.8
% box 9:  0 .8<= PHIT < 0.9
% box 10:  0.9 <= PHIT <= 1 

% This is the total amount of times each individual seen a specified range
% of phit values.
TotalPhitSeen = zeros(50,10);
% phit hit is ranged from 0 to 1 in increments of 0.1
for i = 1:1:50 % iterate over all individuals
    for j = 1:1:16 % iterate over relevant trials
        maxTime = Tm(i,j);
        if isnan(maxTime) ~= 1 % make sure maxTime is a number
            for k = 0:1:maxTime   % iterate over the game (non-evacuated) time of individual i in trial j
                for n = 1:1:10 % Phit bins
                    if n >= 1 && n < 10
                        if Qm(j,k+1) < (n)/10 && Qm(j,k+1) >= (n-1)/10 
                            TotalPhitSeen(i,n)=TotalPhitSeen(i,n)+1; % if the Phit value is within the bin, 
                                                                     % increase the total Phit seen by 1
                        end
                    end
                    if n == 10
                       if Qm(j,k+1) <= (n)/10 && Qm(j,k+1) >= (n-1)/10 
                            TotalPhitSeen(i,n)=TotalPhitSeen(i,n)+1; % if the Phit value is within the bin, 
                                                                     % increase the total Phit seen by 1
                        end 
                    end
                end
            end
        end
    end
end
TotalPhitSeen;
size(TotalPhitSeen); % [50 - by - 10] matrix
length(TotalPhitSeen);

% This is the total amount of times each individial evacuated in a
% specified range of phit values.
TotalPhitEvac = zeros(50,10);
for i = 1:1:50; % iterate over all individials
    for j = 1:1:16; % iterate over all trilas
        for n = 1:1:10; % Phit bins
            if n >= 1 && n < 10
                if Pm(i,j) < (n)/10 && Pm(i,j) >= (n-1)/10
                    TotalPhitEvac(i,n) = TotalPhitEvac(i,n)+1; % if the Phit value is within the bin, 
                                                               % increase the total Phit Evac by 1
                end
            end
            if n == 10
                if Pm(i,j) <= (n)/10 && Pm(i,j) >= (n-1)/10
                    TotalPhitEvac(i,n) = TotalPhitEvac(i,n)+1; % if the Phit value is within the bin, 
                                                               % increase the total Phit Evac by 1
                end
            end
        end
    end
end
TotalPhitEvac;
size(TotalPhitEvac); % [50 - by - 10] matrix
length(TotalPhitEvac);

% This is the probability for each individual to evacuate in the specified
% Phit range {0,0.1,0.2, ... , 0.9,1}
q = TotalPhitSeen.\TotalPhitEvac; % With this q matrix, we will do a best linear fit,
size(q);                          % from that, we will get the average of
length(q);                        % all the slopes then use that value as
                                  % the linear parameter in the linear
                                  % model.

% ############################################################################################################
% ############################################################################################################


% ############################################################################################################
% ############################################################################################################
% Concatinate matrix where we get only a row vector
Cr = [C1(19,:) C1(28,:) C1(29,:) C1(49,:) C1(64,:) C1(67,:) C1(76,:) C1(102,:) C1(108,:)...
    C1(112,:) C1(113,:) C1(123,:) C1(125,:) C1(130,:) C1(143,:) C1(144,:)]; 
Qr = [Q1(19,:) Q1(28,:) Q1(29,:) Q1(49,:) Q1(64,:) Q1(67,:) Q1(76,:) Q1(102,:) Q1(108,:)...
    Q1(112,:) Q1(113,:) Q1(123,:) Q1(125,:) Q1(130,:) Q1(143,:) Q1(144,:)]; 

% we will sort the values of Qr and specify the index
[Qrs, Ir] = sort(Qr);
Crs = Cr(Ir); % this sorts Cr in the same way Qr was sorted

% we will split Qrs into bins of specified lengths then store the Crs
% values in those bins corresponding to the index of Qrs, 
% then make a sparse matrix , take the average of each bin, then make the 
% sparse matrix into a full matrix

% bin the Qrs
numbins1 = 11;
binranges1 = linspace(min(Qrs),max(Qrs),numbins1);
[bincounts1,ind1] = histc(Qrs,binranges1); % bincounts = number of values within each bin 
                                            % ind = bin number each entry 
                                            % of Qrs sorts into
cy = sparse(1:length(Qrs),ind1,Crs); % (Qrs index, Which bin its in) the value of Crs at the Qrs Index
                                      % cy = (Qrs Index, Bin number)  Crs(at Qrs Index)      
a = sum(cy)./sum(cy~=0); % this takes the average of each bin

% this will be our new averaged cumulative evacuate within in the Phit range of 0:0.1:1
Crs_avg = full(a); % This makes the sparse matrix to full storage organization
% ############################################################################################################
% ############################################################################################################



% We wll construct a function which takes the MLE function and maximizes our Hill-function 
% Parameters. (Lambda, K(Threshold Value), n(Steepness of Curve)
% The sample data must be a vector
% phat = mle(data,'distribution',dist) returns parameter estimates for a distribution specified by dist.

% ############################################################################################################
% ############################################################################################################
% we need to construct the total number of individual at home and evacaute
% for each PHIT

% Total number of participant who are ATHOME and observe PHIT
H = zeros(1,10); 

% Total number of participant who have Evacauted and observe PHIT
J = zeros(1,10); 

% Find time which experiment ended "Qm"
TrialEndTime = zeros(1,16); % we will store the time at which the trials ended
for j = 1:16 % iterate through each trial
    for i = 2:61 % iterate through each time step until Qm reaches 0 or 1 
        if Qm(j,i) == 1;
            Qm(j,i);
            TrialEndTime(j) = i; % Time when disaster struck
            j; % Relavent trial to time
            break;
        end
        if Qm(j,i) == 0;
            Qm(j,i);
            TrialEndTime(j) = i; % Time when disaster missed
            j; % Relavent trial to time
            break;
        end
    end
end


for i = 1:1:10 % iterate through the possible PHIT ranges {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1}
    for j = 1:50 % iterate through each individual
        for r = 1:16 % iterate through each trial
            maxTime = TrialEndTime(r); % this is the time the trail (r) ended
            IndvEvacTime = Tm(j,r); % this is the time the individual (j) evacuated in trial (r)
            for k = 1:maxTime % iterate until trial (r) is done
                if k < IndvEvacTime && isnan(IndvEvacTime) ~= 1% iterate until the individual evacuates
                    if i >= 1 && i < 10
                        if Qm(r,k) < i/10 && Qm(r,k) >= (i-1)/10 % this makes sure we are in the correct PHIT range
                            H(i) = H(i) + 1;
                        end
                    end
                    if i == 10
                        if Qm(r,k) <= i/10 && Qm(r,k) >= (i-1)/10 % this makes sure we are in the correct PHIT range
                            H(i) = H(i) + 1;
                        end
                    end
                       
                end
                if k >= IndvEvacTime && k <= maxTime % iterate until indv evacutes and trial (r) ends
                    if i >= 1 && i < 10
                        if Qm(r,k) < i/10 && Qm(r,k) >= (i-1)/10 % this makes sure we are in the correct PHIT range
                            J(i) = J(i) + 1;
                        end
                    end
                    if i == 10
                        if Qm(r,k) <= i/10 && Qm(r,k) >= (i-1)/10 % this makes sure we are in the correct PHIT range
                            J(i) = J(i) + 1;
                        end
                    end
                end
            end
        end
    end
end
H; % this is the number of time for all the trials and all individuals who STAYED AT HOME for PHIT ranges
size(H); % [1 - by - 10]
length(H);

J; % this is the number of time for all the trials and all individuals who EVACUATED for PHIT ranges
size(J); % [1 - by - 10]
length(J);    


% ############################################################################################################
% ############################################################################################################




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
title(['Averaged Cumulative Number Evacauted ', '  N = ', num2str(N-1)])
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
plot(T_range, C)
title(['Cumulative Number Evacauted ', 'Trial ', num2str(z), '  N = ', num2str(N-1)])
ylabel('Number Evacuated');
xlabel('Time (seconds)');
axis([T_range1,0,50]);

% This data value will be compared to the mean number evacuated in the
% theoretical data
%================================================================================================


%================================================================================================
% This is PHIT vs TIME

figure();
plot(T_range, Q);
title(['Phit Trajectory ', 'Trial ', num2str(z), '  N= ',num2str(N-1)]);
ylabel('Probability of Disaster to Strike');
xlabel('Time (seconds)');
axis([T_range1,0,1]);

%================================================================================================

 
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
