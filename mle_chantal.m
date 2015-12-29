% desired trial #s are 19,28,29,49,64,67,76,102,108,112,113,123,125,130,143,144
z = [19,28,29,49,64,67,76,102,108,112,113,123,125,130,143,144];
% misses 28 29 64 67 113 125 130 143 144
% hits 19 49 76 102 108 112 123

% import disaster probabilities (P_hit) 
all_P_hits = importdata('gameinfo.xlsx');

% round P_hits down to nearest tenth, since participants only see discrete
% values of P_hit
round_P_hits = 10.^floor(log10(abs(all_P_hits)));
round_P_hits = fix(all_P_hits./round_P_hits).*round_P_hits;
round_P_hits(isnan(round_P_hits))=0;
round_P_hits = round(round_P_hits,1);

% import number evacuated
all_evac = importdata('evapeocumu.xlsx');

ns = length(z); % number of samples

all_time_at_evac = importdata('evacuateTime.xlsx'); % time at which each individual evacuates
all_Phit_at_evac = importdata('evacuateProb.xlsx'); % Phit at which each individual evacuates

P_hits = zeros(61,length(z));
evac = P_hits;
time_at_evac = zeros(50,length(z));
Phit_at_evac = time_at_evac;

% this is our data set
for i = 1:length(z)
    P_hits(:,i) = round_P_hits(z(i),:);
    evac(:,i) = all_evac(z(i),:); 
    time_at_evac(:,i) = all_time_at_evac(:,z(i));
    Phit_at_evac(:,i) = all_Phit_at_evac(:,z(i));
end

P_hit_range=0:0.1:1;

% Number of times participants at home observe P_hit values
H = zeros(11,length(z)); 

% Number of times participants evacuate at P_hit values
J = zeros(11,length(z)); 

% Find times when game ended 
TrialEndTime = zeros(length(z),1); 
for j = 1:length(z) % iterate through each trial
    for i = 2:61 % iterate through each time step until P_hit is 0 or 1 
        if P_hits(i,j) == 1
            TrialEndTime(j) = i; 
            break;
        end
        if P_hits(i,j) == 0
            TrialEndTime(j) = i; 
            break;
        end
    end
end

for j = 1:50 % iterate through each individual
    for r = 1:length(z) % iterate through each trial
        maxTime = TrialEndTime(r); % this is the time the trial (r) ended
        IndvEvacTime = time_at_evac(j,r); % this is the time the individual (j) evacuated in trial (r)
        for k = 1:maxTime % iterate until trial (r) is done
            if k <= IndvEvacTime && IndvEvacTime ~= -1 % iterate until the individual evacuates in trial (r)
                H(round((P_hits(k,r)+0.1)*10),r)=H(round((P_hits(k,r)+0.1)*10),r)+1;     
            end
            if k == IndvEvacTime % iterate until indv evacuates in trial (r)
                J(round((P_hits(k,r)+0.1)*10),r)=J(round((P_hits(k,r)+0.1)*10),r)+1;     
            end
        end
    end
end

cumuH = sum(H,2);
cumuJ = sum(J,2);

theta = [1; 0.5; 10];
hillfun = @(theta)hill(cumuH(2:end),cumuJ(2:end),theta,P_hit_range(2:end)); % cannot include P_hit = 0
A = [1 0 0; 0 1 0; 0 0 0];
b = [1 1 100];
options = optimoptions(@fminunc,'MaxFunEvals',10000);
[theta_MLE_unc,ML_unc] = fminunc(hillfun,theta,options);
[theta_MLE,~] = fmincon(hillfun,theta,A,b);

