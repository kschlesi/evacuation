% Script for Loss Matrix Analysis.
%% Question: For which loss matrix is the empirical strategy optimal?

% (1) Find best-fit strategy to empirical behavior.
% (2) Compute final scores from actual loss matrix.
% (3) Determine range of test loss matrices to use.
% (4) Compute expected score from each test matrix.
% (5) Plot & compare scores from all loss matrices.

%% Find best-fit strategy to empirical behavior.

% load evac data 
bins = -0.05:0.1:1.05;                  % this gives the bin edges
[~,~,~,~,missing] = load_evac_data(0);  % find missing trials
trials = trial_ix('ind',50,1,missing);  % find all Ind-50 trials
[H,J,Theta,AvgCumEvac,~,Q1,T1,P1,C1] = load_evac_data(0);  % load data
[ntrials,nts] = size(Q1);               % get parameters (trials, timesteps
[N,~] = size(P1);                       % and subjects)

% find best-fit powerlaw strategy
qform = @(Phit_,pv_) pv_(1).*Phit_.^pv_(2);
startp = [1.1;1];
Phit = 0:0.1:1;
params = ML_fit_beta(qform,Phit(2:end-1),H(:,2:end-1),J(:,2:end-1),startp);

%% Compute final scores from actual loss matrix.

[e_score,e_trscore] = total_score(qform,params,Phit,trials,Q1,loss_matrix);

%% Determine range of test loss matrices to use.
%% Compute expected score from each test matrix.

% rho_a = all assets as percentage of span (range: 0:1)
% rho_b = ratio of travel assets to posessions (range: 0:1)
% original loss matrix: he = 6, hn = 10, me = 2, mn = 0
%                       rho_a = 0.6, rho_b = 0.5

rho_a_list = 0:0.1:1;
rho_b_list = 0:0.1:1;
t_scores = zeros(numel(rho_a_list),numel(rho_b_list));
t_trscores = zeros(numel(rho_a_list),numel(rho_b_list),numel(trials));
for rho_a = rho_a_list
  for rho_b = rho_b_list
    disp(['rho_a = ' num2str(rho_a) ', rho_b = ' num2str(rho_b)]);
    [t_score,t_trscore] = total_score(qform,params,Phit,trials,Q1,...
                                      loss_matrix_unitless(rho_a,rho_b));
    t_scores(rho_a_list==rho_a,rho_b_list==rho_b) = t_score;
    t_trscores(rho_a_list==rho_a,rho_b_list==rho_b,:) = shiftdim(t_trscore,-2);
  end
end

%% Compute 

%% Figure

figure; bcolor(t_scores,rho_a_list,rho_b_list); 
        colorbar;
