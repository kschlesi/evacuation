function PbXStochastic_sample_MAIN_2()
% Considering Binomial Distribution
% MODIFICATIONS: We will make q=Phit as a linear function, but Phit is
%                modelled as a stochatic process. Phit values will come
%                from experiment ran values

format long
% this will use a Binomial Distribution to model the transition matrix
% this is still for individual choices not group choices

% this function uses the matlab ODE solver 'ode45' to solve a master
% equation, from a given initial value condition, to find the time
% evolution of the system given a matrix A.

% Time Variable: define range of timesteps in which to give solutions
T_range1 = [0, 60];   % Time evolution is seconds of trial
T_range = 0:1:T_range1(2);
% Parameters: define other parameters in the equation
% Condition on transiiton matix: We assume the time step is extreamly small
%                                that way, only one subject can evauate at
%                                a time and not more then one can evacuate
%                                in that timestep.
% IGNORE CONDITION ON TRANSITION MATRIX SINCE WE CONSIDE  R A DISTRIBUTION
% NOW

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

% Initial Values for State Variables 
P_init = zeros(N,1);     % initial probability of being in states
for n = 1:N 
    if n == 1
        P_init(n) = 1;
    end
    if n > 1
        P_init(n) = 0;
    end
end
size(P_init);
length(P_init);
P_init;

% ******Phit is a random stochastic process******
% make probability of one person to evacuate still linear to Phit but Phit
% is now a stochastic process that is given from the experiment data, so we
% need to import the file probability.csv

% relevent trials are 19,28,29,49,64,67,76,102,108,112,113,123,125,130,143,144
z = 49; % this is the trial number

Q1 = importdata('gameinfo.xlsx'); 
length(Q1);
size(Q1); % row numbers = trial number, column numbers = time (0:1:60)
Q = Q1(z,:); % this is a row vector 
length(Q);  % length(Q1) = total number of time-steps (0 to 60)
size(Q);


% import experiment cumulative number evacauted
C1 = importdata('evapeocumu.xlsx');
length(C1);
size(C1); % row numbers = individual number , column numbers = trial number
C = C1(z,:); % this is a row vector
length(C); % length(C1) = total number of time-steps (0 to 60)
size(C);

% already modified all the data to account for t = 0 and n = 0
% experimental values did not include time t = 0 where Phit = 0.5

% if we need to increase column by 1 and add a column vector of zeros for all trials.
% c = size(Q);
% c1 = c(1); % chooses row length of Q, c(2) would choose column length of Q
% NewCol = zeros(c(1), 1); % new column vector of size [c1 1]
% Q = [NewCol Q]; % new matrix of size row(Q) - by - (col(Q)+1)

% ##############################################################################################################
% ##############################################################################################################

% HILL-FUNCTION MODEL
xH = [104.1,3.6,4.37]; % xH = [lambda,k,n]
qH = xH(1)*Q.^xH(3)./(Q.^xH(3)+xH(2).^xH(3)); 


% LINEAR DECISON MODEL
xL = [0.06014]; % xL = [m (slope)]
qL = xL(1)*Q; % This is a ROW vector 
length(qL);  % length(qL) = total number of time steps
size(qL);

% THRESHOLD MODEL
xT = [0.55,0.8]; % xT = [k,lambda]
qT = zeros(1,length(Q));
for n = 1:1:length(Q)
    if Q(n) >= 0 && Q(n) < xT(1)
        qT(n) = 0;
    end
    if Q(n) >= xT(1) && Q(n) <= 1
        qT(n) = xT(2);
    end
end
qT;
length(qT);  % length(qT) = total number of time steps
size(qT);


% ##############################################################################################################
% ##############################################################################################################


% ##############################################################################################################
% ##############################################################################################################
% make the transition matrix
% let n = rows, let i = columns

% We need A[ n | n ] = - SUM[ i ~= n] P[ i | n ]  ************
% nchoosek(N-i,n-i) = (factorial(N-i)/(factorial(n-i)*factorial(N-n)))


% HILL FUNCTION
AH = zeros(N,N,length(qH)); % make N-by-N matrix of zeros
for m=1:length(qH)
    for n=1:N
      for i=1:N
        if i < n
          AH(n,i,m)=nchoosek(N-i,n-i).*qH(m)^(n-i).*(1-qH(m))^(N-n);
        end
      end
    end
    BH=sum(AH);  % if A is a matrix, then sum(A) returns a row vector containing
               % the sum of each column.
    for i=1:N
        AH(i,i,m)=-BH(1,i,m); % ************
    end
end
size(AH); % N - by - N - by - length(q) matrix
length(AH);
AH;
BH;

% we need to make A into a 2D matrix
% reshape(A, N,N,length(q) = reshape(A, [], length(q)) gives a N*N - by -
% length(q) matrix

AHr = reshape(AH, [], length(qH)).'; % we transpose the shape after by .'
size(AHr); % length(q) - by - N*N
AH = AHr;

% Solver Options
options = [];

% Solver: solve the equation defined in 'odefunc' over the time range,
% starting with initial conditions P_init. this returns the variable T,
% a ts x 1 vector containing a list of times at which solutions are given,
% and P, a ts x length(P_init) vector containing the value of each entry of
% P(t) at each timestep in T.

[T, PH] = ode45(@(t,PH)odefunc(PH, AH, AHr, N, t), T_range, P_init, options);

length(T); % n (rows) " this is from 0 to 60 seconds"
size(T); % n - by - 1 column vector
T;

length(PH); % n (rows) " this is from 0 to 60 seconds"
size(PH);  % Length(T) - by - N  (rows - by - columns)
PH; % matrix of size [length(T) - by - N] ROWS = Time (0 to 60), COLUMNS = Number Evacauted (0 to 50)
% ##############################################################################################################
% ##############################################################################################################


% ##############################################################################################################
% ##############################################################################################################
% make the transition matrix
% let n = rows, let i = columns

% We need A[ n | n ] = - SUM[ i ~= n] P[ i | n ]  ************
% nchoosek(N-i,n-i) = (factorial(N-i)/(factorial(n-i)*factorial(N-n)))


% LINEAR FUNCTION
AL = zeros(N,N,length(qL)); % make N-by-N matrix of zeros
for m=1:length(qL)
    for n=1:N
      for i=1:N
        if i < n
          AL(n,i,m)=nchoosek(N-i,n-i).*qL(m)^(n-i).*(1-qL(m))^(N-n);
        end
      end
    end
    BL=sum(AL);  % if A is a matrix, then sum(A) returns a row vector containing
               % the sum of each column.
    for i=1:N
        AL(i,i,m)=-BL(1,i,m); % ************
    end
end
size(AL); % N - by - N - by - length(q) matrix
length(AL);
AL;
BL;

% we need to make A into a 2D matrix
% reshape(A, N,N,length(q) = reshape(A, [], length(q)) gives a N*N - by -
% length(q) matrix

ALr = reshape(AL, [], length(qL)).'; % we transpose the shape after by .'
size(ALr); % length(q) - by - N*N
AL = ALr;

% Solver Options
options = [];

% Solver: solve the equation defined in 'odefunc' over the time range,
% starting with initial conditions P_init. this returns the variable T,
% a ts x 1 vector containing a list of times at which solutions are given,
% and P, a ts x length(P_init) vector containing the value of each entry of
% P(t) at each timestep in T.

[T, PL] = ode45(@(t,PL)odefunc(PL, AL, ALr, N, t), T_range, P_init, options);

length(T); % n (rows) " this is from 0 to 60 seconds"
size(T); % n - by - 1 column vector
T;

length(PL); % n (rows) " this is from 0 to 60 seconds"
size(PL);  % Length(T) - by - N  (rows - by - columns)
PL; % matrix of size [length(T) - by - N] ROWS = Time (0 to 60), COLUMNS = Number Evacauted (0 to 50)
% ##############################################################################################################
% ##############################################################################################################


% ##############################################################################################################
% ##############################################################################################################
% make the transition matrix
% let n = rows, let i = columns

% We need A[ n | n ] = - SUM[ i ~= n] P[ i | n ]  ************
% nchoosek(N-i,n-i) = (factorial(N-i)/(factorial(n-i)*factorial(N-n)))


% THRESHOLD FUNCTION
AT = zeros(N,N,length(qT)); % make N-by-N matrix of zeros
for m=1:length(qT)
    for n=1:N
      for i=1:N
        if i < n
          AT(n,i,m)=nchoosek(N-i,n-i).*qT(m)^(n-i).*(1-qT(m))^(N-n);
        end
      end
    end
    BT=sum(AT);  % if A is a matrix, then sum(A) returns a row vector containing
               % the sum of each column.
    for i=1:N
        AT(i,i,m)=-BT(1,i,m); % ************
    end
end
size(AT); % N - by - N - by - length(q) matrix
length(AT);
AT;
BT;

% we need to make A into a 2D matrix
% reshape(A, N,N,length(q) = reshape(A, [], length(q)) gives a N*N - by -
% length(q) matrix

ATr = reshape(AT, [], length(qT)).'; % we transpose the shape after by .'
size(ATr); % length(q) - by - N*N
AT = ATr;

% Solver Options
options = [];

% Solver: solve the equation defined in 'odefunc' over the time range,
% starting with initial conditions P_init. this returns the variable T,
% a ts x 1 vector containing a list of times at which solutions are given,
% and P, a ts x length(P_init) vector containing the value of each entry of
% P(t) at each timestep in T.

[T, PT] = ode45(@(t,PT)odefunc(PT, AT, ATr, N, t), T_range, P_init, options);

length(T); % n (rows) " this is from 0 to 60 seconds"
size(T); % n - by - 1 column vector
T;

length(PT); % n (rows) " this is from 0 to 60 seconds"
size(PT);  % Length(T) - by - N  (rows - by - columns)
PT; % matrix of size [length(T) - by - N] ROWS = Time (0 to 60), COLUMNS = Number Evacauted (0 to 50)
% ##############################################################################################################
% ##############################################################################################################


close all
%================================================================================================
%  Calculate Mean Number Evacuated with the P vector
PH1 = PH*((0:(N-1)).'); % size ( Length(T_range) - by - 1)
PL1 = PL*((0:(N-1)).');
PT1 = PT*((0:(N-1)).');

cb = 'kbgrm';
figure();
cla;
plot(T_range, PH1, 'color', cb(1)) % Hill Function BLACK
hold on
plot(T_range, PL1, 'color', cb(2)); % Linear Function BLUE
plot(T_range, PT1, 'color', cb(3)); % Threshold Function GREEN
[hAx,hLine1,hLine2] = plotyy(T_range, C, T_range, Q); % Experiment Number Evacuated RED, Phit Trajectory MAGENTA
% hAx(1) is left y-axis handle, hAx(2) is right y-axis handle, hLine1 is {C} plot, hLine2 is {Q} plot

title([ 'Trial ', num2str(z), '  N = ', num2str(N-1)], 'FontWeight', 'Bold', 'FontSize', 10)
legend('Hill','Linear','Threshold','Experiment','PHIT','Location','NorthWest');
xlabel('Time (seconds)', 'FontWeight', 'Bold', 'FontSize', 10); 
set(get(hAx(1), 'YLabel'), 'String', 'Number Evacuated', 'FontWeight', 'Bold', 'FontSize', 10);% Left y - axis hAx(1)
set(hAx(1), 'YLim', [0,50], 'YTick', [0:5:50], 'YColor' , cb(1)); 
set(get(hAx(2), 'YLabel'), 'String', 'Disaster Probability', 'FontWeight', 'Bold', 'FontSize', 10); % Right y - axis hAx(2)
set(hAx(2), 'YLim', [0,1], 'YTick', [0:0.1:1], 'YColor' , cb(5));
set(hLine1,'color', cb(4), 'LineStyle', '--'); % Cumulative Number Evacuated RED
set(hLine2,'color',cb(5), 'LineStyle', ':'); % Phit Trajectory MAGENTA
%axis([T_range1,0,50,0,1]);

%================================================================================================


%================================================================================================
% gives rows of PH on x-axis and columns of PH on y-axis

% Calculate the standard deviation at each time step
% At each timestep, subtract the mean value from the vector 0:50, then square it, 
% multiply by the probability distribution at the timestep or divide by (N-1), and take the square root.

% Calculate Standard Deviation
% make number evacuted vector
V = 0:1:N-1; % column vector of size ( 1 - by - N)
sumH = 0;
sumH1 = zeros(length(T_range),1);
stdevH = zeros(length(T_range),1); % standard deviation for every time step, size (length(T_range) - by - 1)
PHplus = zeros(length(T_range),1); % 2 standard deviations from mean
PHminus = zeros(length(T_range),1); % 2 standard deviations from mean

for i = 1:1:length(T_range) % Iterate through all time steps
    for j = 1:1:N % Iterate through all number of individuals evacauted
        sumH = V(j) - PH1(i);
        sumH1(i) = sumH^2;
        %stdevH(i) = sqrt(PH(i)*sumH1(i)); 
        stdevH(i) = sqrt(sumH1(i)/(N-1));
    end 
    PHplus(i) = PH1(i) + 2*stdevH(i);
    PHminus(i) = PH1(i) - 2*stdevH(i);
end

ca = 'kbgrm';
figure('Position',[1 1 1000 800]);
bcolor(PH');
hH = colorbar;
hold on
plot(C, '--', 'color', ca(5));
plot(T_range, PHplus, ':', 'color', 'y');
plot(T_range, PHminus, ':', 'color', 'y');

title(['Hill ', ' Trial ', num2str(z), '  N = ', num2str(N-1)], 'FontWeight', 'Bold', 'FontSize', 10);
legend('Hill','Experiment','Location','NorthWest');
ylabel('Number Evacuated', 'FontWeight', 'Bold', 'FontSize', 10);
xlabel('Time (seconds)', 'FontWeight', 'Bold', 'FontSize', 10);
set(get(hH, 'Title'), 'String', 'Evacuation Probability', 'FontWeight', 'Bold', 'FontSize', 10,'Position',[700 700]); 
LabelPosH = get(get(hH,'Title'),'Position');
set(get(hH,'Title'),'Position',[LabelPosH(1) + 6.2 LabelPosH(2) - 0.5  LabelPosH(3)],'Rotation',90);


% gives rows of PL on x-axis and columns of PL on y-axis

% Calculate Standard Deviation
% make number evacuted vector
V = 0:1:N-1; % column vector of size ( 1 - by - N)
sumL = 0;
sumL1 = zeros(length(T_range),1);
stdevL = zeros(length(T_range),1); % standard deviation for every time step, size (length(T_range) - by - 1)
PLplus = zeros(length(T_range),1); % 2 standard deviations from mean
PLminus = zeros(length(T_range),1); % 2 standard deviations from mean

for i = 1:1:length(T_range) % Iterate through all time steps
    for j = 1:1:N % Iterate through all number of individuals evacauted
        sumL = V(j) - PL1(i);
        sumL1(i) = sumL^2;
        %stdevL(i) = sqrt(PL(i)*sumL1(i)); 
        stdevL(i) = sqrt(sumL1(i)/(N-1));
    end 
    PLplus(i) = PL1(i) + 2*stdevL(i);
    PLminus(i) = PL1(i) - 2*stdevL(i);
end


function h = bcolor(inmat)
% provides a balanced color plot (no row/cols left out) with no edge lines
    if ~ismatrix(inmat)
        error('input matrix must be two-dimensional'); 
    end
    
    pad = mean(mean(inmat));
    h = pcolor(padarray(inmat,[1 1],pad,'post'));
    set(h, 'EdgeColor', 'none');
    
end


ca = 'kbgrm';
figure();
bcolor(PL');
hL = colorbar;
hold on
plot(C, '--', 'color', ca(5));
plot(T_range, PLplus, ':', 'color', 'y');
plot(T_range, PLminus, ':', 'color', 'y');

title(['Linear ', ' Trial ', num2str(z), '  N = ', num2str(N-1)], 'FontWeight', 'Bold', 'FontSize', 10);
legend('LInear','Experiment','Location','NorthWest');
ylabel('Number Evacuated', 'FontWeight', 'Bold', 'FontSize', 10);
xlabel('Time (seconds)', 'FontWeight', 'Bold', 'FontSize', 10);
set(get(hL, 'Title'), 'String', 'Evacuation Probability', 'FontWeight', 'Bold', 'FontSize', 10); 
LabelPosL = get(get(hL,'Title'),'Position');
set(get(hL,'Title'),'Position',[LabelPosL(1) + 6.2 LabelPosL(2) - 0.5  LabelPosL(3)],'Rotation',90);

% gives rows of PT on x-axis and columns of PT on y-axis

% Calculate Standard Deviation
% make number evacuted vector
V = 0:1:N-1; % column vector of size ( 1 - by - N)
sumT = 0;
sumT1 = zeros(length(T_range),1);
stdevT = zeros(length(T_range),1); % standard deviation for every time step, size (length(T_range) - by - 1)
PTplus = zeros(length(T_range),1); % 2 standard deviations from mean
PTminus = zeros(length(T_range),1); % 2 standard deviations from mean

for i = 1:1:length(T_range) % Iterate through all time steps
    for j = 1:1:N % Iterate through all number of individuals evacauted
        sumT = V(j) - PT1(i);
        sumT1(i) = sumT^2;
        %stdevT(i) = sqrt(PT(i)*sumT1(i)); 
        stdevT(i) = sqrt(sumT1(i)/(N-1));
    end 
    PTplus(i) = PT1(i) + 2*stdevT(i);
    PTminus(i) = PT1(i) - 2*stdevT(i);
end

ca = 'kbgrm';
figure();
bcolor(PT');
hT = colorbar;
hold on
plot(C, '--', 'color', ca(5));
plot(T_range, PTplus, ':', 'color', 'y');
plot(T_range, PTminus, ':', 'color', 'y');

title(['Threshold ', ' Trial ', num2str(z), '  N = ', num2str(N-1)],'FontWeight', 'Bold', 'FontSize', 10);
legend('Threshold','Experiment','Location','NorthWest');
ylabel('Number Evacuated', 'FontWeight', 'Bold', 'FontSize', 10);
xlabel('Time (seconds)', 'FontWeight', 'Bold', 'FontSize', 10);
set(get(hT, 'Title'), 'String', 'Evacuation Probability', 'FontWeight', 'Bold', 'FontSize', 10); 
LabelPosT = get(get(hT,'Title'),'Position');
set(get(hT,'Title'),'Position',[LabelPosT(1) + 6.2 LabelPosT(2) - 0.5  LabelPosT(3)],'Rotation',90);


%================================================================================================


end  % end main function


function dP = odefunc(P, A, Ar, N, t)
% This function defines the differential equation that will be called in
% the main function. The first two arguments MUST be t (the timestep 
% calculated internally by ode45) and the variable being solved for (P).
% Other parameters (A,qt,Ar,N,t) must be passed as arguments.
% The output should be the equivalent of dP/dt for one timestep.

Ari = interp1(0:(size(Ar,1)-1),Ar,t); % interpolate data to get value at specified time
A = reshape(Ari,N,N); % returns a N - by - N matrix whose elements are taken 
                      % column - wise from Ari
dP = A*P;  

end