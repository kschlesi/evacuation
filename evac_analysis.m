%%% script for running things with evacuation data
addpath(genpath(pwd));

% load_evac_data 
bins = 0:0.1:1.1;
%bins = -0.05:0.1:1.05;  % this gives the bin EDGES
trials = [19,28,29,49,64,67,76,102,108,112,113,123,125,130,143,144];
[H,J,Theta,AvgCumEvac,missing,Q1,T1,P1,C1] = load_evac_data(0,trials,bins);
[ntrials,n_ts] = size(Q1);
[N,~] = size(P1);

%% we can (a) fit experimental data to a particular q function 
% with maximum likelihood estimation

Phit = 0:0.1:1;
% qform = @(Phit_,pv_) pv_(1).*(Phit_.^pv_(3))./(pv_(2).^pv_(3)+Phit_.^pv_(3));
% startp = [1;0.5;10];
qform = @(Phit_,pv_) pv_(1).*Phit_.^pv_(2);
%startp = [1.1,1];
startp = [0.8;10];
params = ML_fit_beta(qform,Phit(:,2:end),H(:,2:end),J(:,2:end),startp);

figure(1); hold all; plot(Phit,qform(Phit,params));
xlabel('Phit'); ylabel('evacuation probability'); title('best-fit hill model, ind. trials, shelterSpace=50');
%legend(['\Lambda = ' num2str(params(1)) ', k = ' num2str(params(2)) ', n = ' num2str(params(3))]);
%legend(['\alpha = ' num2str(params(1)) ', \beta = ' num2str(params(2))]);
hold off;

%% fit for all group games ss=50

% trials = unique([trial_ix('fTG',50,5,missing);...
%             trial_ix('lTG',50,5,missing);...
%             trial_ix('mR',50,5,missing)])';
% trials = unique([trial_ix('fTG',50,25,missing);...
%            trial_ix('lTG',50,25,missing);...
%            trial_ix('mR',50,25,missing)])';
% trials = unique([trial_ix('fTG',50,5,missing);...
%             trial_ix('fTG',50,25,missing)])';
% trials = unique([trial_ix('lTG',50,5,missing);...
%             trial_ix('lTG',50,25,missing)])';
trials = unique([trial_ix('mR',50,5,missing);...
            trial_ix('mR',50,25,missing)])';
        
[H,J,~,~,~,Q1] = load_evac_data(0,trials,bins);
startp = params;
params = ML_fit_beta(qform,Phit(:,2:end-1),H(:,2:end-1),J(:,2:end-1),startp);

figure(1); hold all; plot(Phit,qform(Phit,params),'--');
xlabel('Phit'); ylabel('evacuation probability'); title('best-fit hill model, ind. trials, shelterSpace=50');
%legend(['\Lambda = ' num2str(params(1)) ', k = ' num2str(params(2)) ', n = ' num2str(params(3))]);
%legend(['\alpha = ' num2str(params(1)) ', \beta = ' num2str(params(2))]);
hold off;        

%% we can (b) solve the master equation for a given trial and q function
%for tr=trials
tr=38;
%L = params(1);
%k = params(2);
%n = params(3);
Phit = 0:0.1:1;
q = @(Phit_) qform(Phit_,params);
[T, P] = solve_master_binom(N,q,Q1(tr,:),Phit);
figure; subplot(1,3,1);
        bcolor(P',T,0:1:N); colorbar; hold on;
        plot(C1(tr,:));
        xlabel('time'); ylabel('# evacuated'); title('cumulative evacuations');
        subplot(1,3,2);
        plot(0:1:59,Q1(tr,:)); axis([0 59 0 1]);
        xlabel('time'); ylabel('Phit'); title('disaster trajectory');
        subplot(1,3,3);
        plot(0:0.1:1,q(0:0.1:1)); axis([0 1 0 1]);
        xlabel('Phit'); ylabel('prob. of each participant to evacuate'); title('decision model');        
%end

%% we can (c) sample trajectories given a particular q-function
% a = vector of propensity functions
% we have ... 1275 possible reactions
% first make a FUNCTION (a) which returns a propensity MATRIX
for tr=trials
%tr = 29; 
nS = 10;
%L = params(1);
%k = params(2);
%n = params(3);
q = @(Phit_) qform(Phit_,params);
A = makeA(N,q(Phit));
tbins = 0:1:60;
% for below, theta_ is a cell array s.t. theta_{1} = A_ and theta_{2} = Phit_traj_
a = @(t_,P_,theta_) a_vec(t_,P_,theta_{1},theta_{2}); % this returns a prop. vector of only lower-triangular vals
aix = tril(ones(N+1),-1);  % indices of lower-triangular entries
theta = {A; Q1(tr,:)}; % these are the non-explicitly-time-dependent input parameters
P0 = [1;zeros(N,1)];   % this is the initial condition
t0 = 0;                % starting time of simulation
endx = find(Q1(tr,:)==Q1(tr,end),1,'first');
T = endx + 1 - (endx==length(Q1(tr,:)));  % duration of simulation
% this creates the update matrix, nu...
% for entry (i,j) -- that means a transition from state idx j to state idx i.
% thus, nu(j,[i,j]) = -1 and nu(i,[i,j]) = +1; convert to vector form (51x1275)
[ix,jx] = ind2sub(size(aix),find(~~aix));
sx = sub2ind(size(aix),find(~~aix));
nu = sparse(jx,1:length(sx),-1,N+1,length(sx)) + sparse(ix,1:length(sx),1);
% this generates the samples with Gillespie algorithm...
[t_mat,n_mat,r_mat,ns_mat,br_mat] = gillespie_tdep(theta,a,nu,P0,t0,T,tbins,nS);

% plot cumulative evacuated
    figure; plot(C1(tr,:),'k--'); hold on; % observed in experiment
            plot(Q1(tr,:).*N,':');         % Phit trajectory
            for i=1:nS
                plot(t_mat(r_mat(:,i)~=0,i),ix(r_mat(r_mat(:,i)~=0,i))-1,'-o'); hold all;
            end
            title(['trial ' num2str(tr) ', ' num2str(nS) ' samples']); axis([0 60 0 50])
            legend('empirical data','Phit trajectory (scaled by N)','location','northwest');
            xlabel('time'); ylabel(['cumulative no. evacuated']);

end  % end loop over trials

%% OR, we can do samples of individuals' trajectories...
% NOTE this needs to be rewritten to take into account discretized Phits
% ^^^ DONE (SEE gillespie_tdep.m)
tr = 112;
nS_ind = 100;
%L = params(1);
%k = params(2);
%n = params(3);
tbins = 0:1:60;
P0_ind = 0;            % initial value of "evacuated" state
t0 = 0; T = 60;        % starting time and duration of simulation
q_ind = @(Phit_) qform(Phit_,params);
% propensity for individual evacuation
a_ind = @(t_,P_,theta_ind_) theta_ind_{1}(floor(theta_ind_{2}(floor(t_)+1)*10)/10)*(1-P_); % there is only one reaction: evacuation
% theta_ind is a cell array such that entry {1} is the fcn q_ind(Phit_)
% and entry {2} is the Phit trajectory for the desired trial to simulate
theta_ind = {q_ind,Q1(tr,:)};
% update matrix: 1 state (evacuated, starts at P_=0), 1 reaction (P_=0 -> 1)
nu_ind = 1;
% call simulator nS times; returns time of evac for each individual
[t_mat,~,r_mat] = gillespie_tdep(theta_ind,a_ind,nu_ind,P0_ind,t0,T,tbins,nS_ind);

% plot hist of evac times & Phit values
    figure; subplot(1,2,1); hist(t_mat(~~r_mat)); hold on;
            ylabel('number of evacuations'); xlabel('evacuation time');
            Phit_traj = [Q1(tr,1),Q1(tr,:)];
            subplot(1,2,2); hist(Phit_traj(floor(t_mat(~~r_mat))+1));
            ylabel('number of evacuations'); xlabel('Phit value');            
            suptitle(['evacuations for ' num2str(nS) ' individuals, trial ' num2str(tr)]);
    
            
%% solve individual amster equation...

% for each trial, solve master equation for 1 individual (2 states)
tr = 1;
%L = params(1);
%k = params(2);
%n = params(3);
Phit = 0:0.1:1;
figure; 
for beta_ = [6,6.2,6.4]
params = [1;beta_];
q = @(Phit_) qform(Phit_,params);
[T, P] = solve_master_binom(1,q,Q1(tr,:),Phit);
        plot(T,P(:,end)'); hold on; 
end
        plot(1:1:60,Q1(tr,:),'k--');
        legend('\beta = 6.0','\beta = 6.2','\beta = 6.4','Phit','location','northwest');
        xlabel('time'); ylabel('prob. of evacuation');
        title(['Individual Evacuation Probability, Trial ' num2str(tr)]);

figure; 
for beta_ = [6,6.2,6.4]
    plot(Phit,qform(Phit,[1;beta_])); hold on;
end
    legend('\alpha = 1, \beta = 6.0','\alpha = 1, \beta = 6.2','\alpha = 1, \beta = 6.4');
    xlabel('Phit'); ylabel('ind. evacuation probability');
    title('Potential Strategies, power-law model');
        
%% fit for STATIC OPTIMAL behavior
% for a given set of trials ('trials')
[Qu,trials] = unique_rows(Q1);
lossmat = loss_matrix(6,10,2,0);
%qform = @(Phit_,pv_) pv_(1).*(Phit_.^pv_(3))./(pv_(2).^pv_(3)+Phit_.^pv_(3));
qform = @(Phit_,pv_) pv_(1).*Phit_.^pv_(2);
Phit = 0:0.1:1; %Phit = Phit(2:end);
% for a single trial: as an example
ptry = params;
scr = exp_score(qform,ptry,Phit,Q1(tr,:),lossmat);

% for all trials:
[tscr,scrs] = total_score(qform,ptry,Phit,trials,Q1,lossmat);

% doing the fit:
%startp = params;
%startp = sopt_params_con;
startp = index_each_cell(prange,optix);
assert(size(startp,1)==1);
SFit = @(ptry_) -1*total_score(qform,ptry_,Phit,trials,Q1,lossmat);
%options = optimoptions(@fminunc,'MaxFunEvals',10000);
%[sopt_params,sopt_score] = fminunc(SFit,startp,options);
A = zeros(length(startp));
b = zeros(length(startp),1);
% constrains all prob. as a function of Phit to be <=1
q_con = @(pvec_)(qform(Phit,pvec_)-1);
[sopt_params_con,sopt_score_con] = fmincon(SFit,startp,A,b,A,b,0,100,@(pvec_)q_add_eq(pvec_,q_con));
%params = sopt_params;

% ok, mapping over the range:
prange = cell(numel(params),1);
prange{1} = linspace(0.6,2,5);
prange{2} = linspace(0.6,1.1,5);
prange{3} = linspace(5,35,8);
psize = cellfun(@numel,prange)';
scoremat = zeros(psize);
tic;
for px = 1:numel(scoremat)
      pix = ind2sub_var_dim(size(scoremat),px);
      par = index_each_cell(prange,pix);
      if any(q_con(par)>0)
          scoremat(px) = NaN;
      else
          scoremat(px) = SFit(par);
      end
end
toc; % about 47 sec for 200 points
best = min(scoremat(:));
allscores = unique(scoremat);
%%%%%%%%%%%%%%% stop here
optix = ind2sub_var_dim(size(scoremat),find(scoremat==best));
optix2 = ind2sub_var_dim(size(scoremat),find((scoremat>allscores(1)).*(scoremat<=allscores(2))));
optix3 = ind2sub_var_dim(size(scoremat),find((scoremat>=allscores(15)).*(scoremat<allscores(35))));
optix4 = ind2sub_var_dim(size(scoremat),find((scoremat>=allscores(35)).*(scoremat<allscores(55))));
optix5 = ind2sub_var_dim(size(scoremat),find((scoremat>=allscores(55)).*(scoremat<allscores(75))));
% plot figures of all scores and max scores
% index_each_cell(prange,optix2)
% figure; scatter3(index_each_cell(prange,optix))%,'r'); hold on;
%         scatter3(index_each_cell(prange,optix2),'b');
%         scatter3(Lrange(Lix3),krange(kix3),nrange(nix3),'k');
%         scatter3(Lrange(Lix4),krange(kix4),nrange(nix4),'g');
%         scatter3(Lrange(Lix5),krange(kix5),nrange(nix5),'y');
%         title('best score parameters');
%         xlabel('L'); ylabel('k'); zlabel('n');
%         set(gca,'XTick',Lrange,'YTick',krange,'ZTick',nrange);
%         axis([Lrange(1),Lrange(end),krange(1),krange(end),nrange(1),nrange(end)])
%         legend([num2str(best) ', max score'],num2str(allscores(15)),...
%                             num2str(allscores(35)),num2str(allscores(55)),...
%                             num2str(allscores(75)),'location','northeast');
        
% what were the mistakes on optimal scores?
all_opt = zeros(size(optix,1),numel(trials));
for i=1:size(optix,1)
    sumall = scoremat(optix(i,:));
    ptry = index_each_cell(prange,optix(i,:));
    [~,all] = total_score(qform,ptry,Phit,trials,Q1,lossmat);
    all_opt(i,:) = all;
end
figure; bcolor([Q1(trials,end)';all_opt]); colorbar;

% plot maximum score found by mesh search of param space
figure;
leg = cell(size(optix,1),1);
for i=1:size(optix,1)
%for i = find(krange(bx(:,2))~=0)
    plot(Phit,qform(Phit,index_each_cell(prange,optix(i,:)))); hold all; 
    %leg{find(krange(bx(:,2))~=0)==i} = ['\Lambda = ' num2str(Lrange(bx(i,1))) ', k = ' num2str(krange(bx(i,2))) ', n = ' num2str(nrange(bx(i,3)))];
    leg{i} = ['\Lambda = ' num2str(prange{1}(optix(i,1))) ', k = ' num2str(prange{2}(optix(i,2))) ', n = ' num2str(prange{3}(optix(i,3)))];
end
xlabel('Phit'); ylabel('Evacuation Probability'); legend(leg,'location','northwest');
title(['Static Optimal Evacuation Strategy (Hill Model), <score> = ' num2str(best)]);

% plot maximum score found overall
figure;
    plot(Phit,qform(Phit,[sopt_params_con(1),sopt_params_con(2),sopt_params_con(3)])); hold all; 
    leg = ['\Lambda = ' num2str(sopt_params_con(1)) ', k = ' num2str(sopt_params_con(2)) ', n = ' num2str(sopt_params_con(3))];
xlabel('Phit'); ylabel('Evacuation Probability'); legend(leg,'location','northwest');
title(['Static Optimal Evacuation Strategy (Hill Model), <score> = ' num2str(sopt_score_con)]);

%% fit for BAYESIAN OPTIMAL behavior
% this updates in the order that trials were presented. use ALL trials.
trials = 1:1:size(Q1,1); trials = trials(:)';
%trials = [19,28,29,49,64,67,76,102,108,112,113,123,125,130,143,144];
                       
% decision model, fcn of Phit and parameters to fit
qhill = @(Phit_,L_,k_,n_) L_.*(Phit_.^n_)./(k_.^n_+Phit_.^n_);
Phit = 0:0.1:1; % list of Phits at which q will be evaluated
q_con = @(pvec_)(qhill(Phit,pvec_(1),pvec_(2),pvec_(3))-1); % constraint: q_con<=0

% set search area of parameter space
Lrange = linspace(0.6,2,5); nL = numel(Lrange);
krange = linspace(0.6,1.1,5); nk = numel(krange);
nrange = linspace(5,35,8); nn = numel(nrange);

% initialize matrices (normalizing 'prob' happens later)
prob = zeros(nL,nk,nn,numel(trials)+1); % posterior probs at each timestep
maxprob = zeros(numel(trials)+1,1);     % maximum normalized posterior prob
n_max_strategy = zeros(numel(trials)+1,1); % number of optimal strategies
prob(:,:,:,1) = ones(nL,nk,nn);         % set flat prior (unnormalized)
n_con = 0;                              % track number of constrained pts
for L=Lrange
  for k=krange
    for n=nrange 
    if any(q_con([L,k,n])>0) % check each proposed strategy with constraint
      prob(Lrange==L,krange==k,nrange==n,1) = 0;
      n_con = n_con+1;
    end
    end
  end
end
maxprob(1) = 1/(nn*nk*nL-n_con);        % set initial normalized max prob.
n_max_strategy(1) = nn*nk*nL-n_con;     % init: all valid strategies optimal

% iterate through each trial and update posterior prob of max score
for tr = trials
disp(['trial ' num2str(tr) ', ' num2str(find(trials==tr)) ' of ' num2str(numel(trials))]);
didHit = Q1(tr,end);
ss = N;
tic;
L_maxscore = zeros(nL,nk,nn);           % likelihood of achieving max score
for L=Lrange
  for k=krange
    for n=nrange
          [~,pEvac] = exp_score(qhill,[L,k,n],Phit,Q1(tr,:),loss_matrix);  
          if didHit % likelihood = pEvac
              L_maxscore(Lrange==L,krange==k,nrange==n) = pEvac;
          else      % likelihood = (1-pEvac)
              L_maxscore(Lrange==L,krange==k,nrange==n) = 1-pEvac;
          end
    end
  end
end
toc; % about 3 sec for 200 parameter space points (one trial)
% update posteriors, normalized max posterior, and # optimal strategies
prob(:,:,:,find(trials==tr)+1) = prob(:,:,:,trials==tr).*L_maxscore;
maxprob(find(trials==tr)+1) = max(unique(prob(:,:,:,find(trials==tr)+1)))./...
                            sum(sum(sum(prob(:,:,:,find(trials==tr)+1))));
n_max_strategy(find(trials==tr)+1) = numel(find(prob(:,:,:,find(trials==tr)+1)==max(unique(prob(:,:,:,find(trials==tr)+1)))));
end
% find final optimal strategy
[Lox,kox,nox] = ind2sub(size(prob(:,:,:,end)),find(prob(:,:,:,end)==max(unique(prob(:,:,:,end)))));
Lopt = Lrange(Lox); kopt = krange(kox); nopt = nrange(nox);

% find strategy evolution
Lt = zeros(sum(n_max_strategy==1),1); kt = Lt; nt = Lt; trt = Lt;
for tr = 1:sum(n_max_strategy==1)-1
  [Loxi,koxi,noxi] = ind2sub(size(prob(:,:,:,find(trials==tr)+1)),find(prob(:,:,:,find(trials==tr)+1)==max(unique(prob(:,:,:,find(trials==tr)+1)))));
  try assert(numel(Loxi)==1)
  catch
      continue
  end
  Lt(tr) = Lrange(Loxi); kt(tr) = krange(koxi); nt(tr) = nrange(noxi); trt(tr) = trials(tr);
end
Lt(end) = Lopt; kt(end) = kopt; nt(end) = nopt; trt(tr) = trials(end);

% plot colorplot of strategy evolution
strategies = zeros(numel(trials),numel(Phit));
for tr = trials
    strategies(trials==tr,:) = qhill(Phit,Lt(trials==tr),kt(trials==tr),nt(trials==tr));
end
figure; bcolor(strategies,Phit,trials); c = colorbar; title('Bayesian strategy evolution');
ylabel('trial'); xlabel('Phit value'); ylabel(c,'evacuation probability');

% plot progression of optimal strategy
trplot = [1;find(~~(diff(Lt)+diff(kt)+diff(nt)))+1];
leg = cell(numel(trplot)+1,1);
leg{1} = ['Final: \Lambda = ' num2str(Lopt) ', k = ' num2str(kopt) ', n = ' num2str(nopt)];
figure; plot(Phit,qhill(Phit,Lopt,kopt,nopt),'-k'); hold all;
for tr=trplot'
    colorval = 1-find(all(bsxfun(@eq,unique_rows(strategies),qhill(Phit,Lt(tr),kt(tr),nt(tr)))'))/size(unique_rows(strategies),1);
    plot(Phit,qhill(Phit,Lt(tr),kt(tr),nt(tr)),'--','Color',[0,colorval,colorval]);
    leg{find(trplot==tr)+1} = ['Trial ' num2str(tr) ': \Lambda = ' num2str(Lt(tr)) ', k = ' num2str(kt(tr)) ', n = ' num2str(nt(tr))];
end
xlabel('Phit'); ylabel('evacuation prob'); title('Bayesian strategy evolution');
legend(leg,'location','northwest'); 

chg = find(sum(diff(strategies),2))+1;
for i=chg'
    figure(1); hold on;
    colorval = find(chg==i)/numel(chg);
    plot(Q1(i-1,:),'--','Color',[0,colorval,colorval]);
    title('trials preceding changes to Bayesian optimal strategy');
    hold off;
end
for i=1:48
    figure(2); hold on;
    colorval = i/48;
    plot(Q1(i,:),'--','Color',[0,colorval,colorval]);
    title('all unique trials');
    hold off;
end

%%%% to do still: make plots of highest probability areas over time

%% HOW does the Bayesian player score? 
lossmat = @(didHit,didEvac) 10*(didHit.*~didEvac) + ...
                             6*(didHit.*didEvac) + ...
                             0*(~didHit.*~didEvac) + ...
                             2*(~didHit.*didEvac) ;
%qform = @(Phit_,pv_) pv_(1).*(Phit_.^pv_(3))./(pv_(2).^pv_(3)+Phit_.^pv_(3));
qform = @(Phit_,pv_) pv_(1).*Phit_.^pv_(2);
Phit = 0:0.1:1; 
q_con = @(pvec_)(qform(Phit,pvec_)-1); % constraint: q_con<=0

%%%%individual Bayesian player:
all_scr = zeros(numel(trials),1);
% choose random strategy to start
ptry = [Lrange(randi(length(Lrange))),krange(randi(length(krange))),nrange(randi(length(nrange)))];
while any(q_con(ptry)>0) % reject if strategy gives prob > 1
    ptry = [Lrange(randi(length(Lrange))),krange(randi(length(krange))),nrange(randi(length(nrange)))];
end
for tr = trials
    % play trial
    all_scr(trials==tr) = exp_score(qform,ptry,Phit,Q1(tr,:),lossmat);
    % update strategy
    ptry = [Lt(trials==tr),kt(trials==tr),nt(trials==tr)];
end
figure; leg = cell(numel(find(all_scr<-9)),1);
for tr = find(all_scr<-9)'
    plot(Q1(tr,:)+0.0005*tr); hold on; 
    leg{find(all_scr<-9)==tr} = ['trial ' num2str(tr)];
end
legend(leg,'location','northwest'); xlabel('time'); ylabel('Phit');
title('No-Evacuate Disaster Hits for Bayesian Optimal Player');
%% COMPARISONS to optimal behavior : part (1), cumulative evacuations
%%%%%%%%%%%%%% and part (2), individual probs and part (3), decision models
%%%% ENSURE CORRECT PARAMS VALUES AND BAYESIAN STRATEGY MATRIX %%%%%%%%%%%%
%trials = 1:48;
trials = trial_ix('ind',50,1,missing);
plot_trials = removeval(trials,missing);
%%%% all Bayesian players (as trained above on ALL trials):
% plotted and compared to expected static optimal, observed cumulative evac

% solve 'static optimal' for comparison
load('score_min_tryALL_cons1.mat');
q_sopt = @(Phit_) qform(Phit_,[sopt_params_con(1),sopt_params_con(2),sopt_params_con(3)]);

% choose random starting strategy from flat prior
ptry = [Lrange(randi(length(Lrange))),krange(randi(length(krange))),nrange(randi(length(nrange)))];
while any(q_con(ptry)>0) % reject if strategy gives prob > 1
    ptry = [Lrange(randi(length(Lrange))),krange(randi(length(krange))),nrange(randi(length(nrange)))];
end
for tr=trials'
    qb = @(Phit_) qform(Phit_,ptry);
    [T, P] = solve_master_binom(N,qb,Q1(tr,:),Phit);
    Cexp = sum(bsxfun(@times,P,0:1:N),2);
    q = @(Phit_) qform(Phit_,params);
    [Tfit,Pfit] = solve_master_binom(N,q,Q1(tr,:),Phit);
    [Tso,Pso] = solve_master_binom(N,q_sopt,Q1(tr,:),Phit);
    C1exp = sum(bsxfun(@times,Pso,0:1:N),2);
    if ismember(tr,plot_trials) && ~ismember(tr,missing)
    figure; subplot(1,2,1);
                bcolor(Pfit',Tfit,0:1:N); colorbar; hold on;
                plot(T,Cexp,'r');
                plot(Tso,C1exp,'y');
                plot(C1(tr,:),'g');
                plot(Q1(tr,:).*N,'--w'); 
                xlabel('time'); ylabel('# evacuated'); title('cumulative evacuations');
                legend('Best-fit','Expected Bayesian Optimal','Expected Static Optimal',...
                       'Empirical','Disaster Trajectory','location','northwest');
%             subplot(1,3,2);
%                 plot(0:1:59,Q1(tr,:),'--k'); axis([0 59 0 1]); hold on;
%                 [Tind, Pind] = solve_master_binom(1,q,Q1(tr,:),Phit);
%                 plot(Tind,Pind(:,end)');
%                 [Tind, Pind] = solve_master_binom(1,qb,Q1(tr,:),Phit);
%                 plot(Tind,Pind(:,end)');
%                 [Tind, Pind] = solve_master_binom(1,q_sopt,Q1(tr,:),Phit);
%                 plot(Tind,Pind(:,end)');
%                 legend('disaster trajectory','best-fit','Bayesian optimal','static optimal','location','northwest');
%                 xlabel('time'); ylabel('probability'); title('individual evacution probabilities');
            subplot(1,2,2);
                plot(0:0.1:1,q(0:0.1:1)); axis([0 1 0 1]); hold on;
                plot(0:0.1:1,qb(0:0.1:1));
                plot(0:0.1:1,q_sopt(0:0.1:1));
                legend('best-fit','Bayesian optimal','static optimal','location','northwest');
                xlabel('Phit'); ylabel('prob. of each participant to evacuate'); title('decision model');        
      set(gcf,'paperunits','inches')
      set(gcf,'papersize',[18,6])
      set(gcf,'paperposition',[0.25,0.25,17.5,5.5])
      %saveas(gcf,['ind_optimal_2pane_tr' num2str(tr) '.jpg']);
    end
    ptry = [Lt(trials==tr),kt(trials==tr),nt(trials==tr)];
end

%% what were actual scores?
trials = 1:48;
lossmat = @(didHit,didEvac) 10*(didHit.*~didEvac) + ...
                             6*(didHit.*didEvac) + ...
                             0*(~didHit.*~didEvac) + ...
                             2*(~didHit.*didEvac) ;
didHits = repmat(Q1(:,end)',N,1);
didEvacs = ~isnan(T1);
scores = sum(lossmat(didHits(:,trials),didEvacs(:,trials)),2);

%% GROUP PREDICTION
shelterSpace = 50;
groupSize = 5;
groupProtocol = 'mR';
p = 100;  % number of times to simulate

% train model on individual trials (output 'params')
trials = trial_ix('ind',shelterSpace,1,missing);
bins = -0.05:0.1:1.05; [H,J,~,~,~,Q1] = load_evac_data(0,trials,bins); n_ts = size(Q1,2);
% qform = @(Phit_,pv_) pv_(1).*(Phit_.^pv_(3))./(pv_(2).^pv_(3)+Phit_.^pv_(3));
% Phit = 0:0.1:1;
% qfit = @(pv_) qform(Phit(:,2:end),pv_);
% MLfit = @(pvec) -1*sum(((sum(H(:,2:end))-sum(J(:,2:end))).*log(1-qfit(pvec)) + sum(J(:,2:end)).*log(qfit(pvec))));
% startp = [1;0.5;10];
qform = @(Phit_,pv_) pv_(1).*Phit_.^pv_(2);
Phit = 0:0.1:1;
qfit = @(pv_) qform(Phit(:,2:end-1),pv_);
MLfit = @(pvec_) -1*sum(((sum(H(:,2:end-1))-sum(J(:,2:end-1))).*log(1-qfit(pvec_)) + sum(J(:,2:end-1)).*log(qfit(pvec_))));
startp = [1.1;1];
options = optimoptions(@fminunc,'MaxFunEvals',10000);
[params1,MLval] = fminunc(MLfit,startp,options);
A = zeros(length(startp)); b = zeros(1,length(startp)); q_con = @(pvec_)(qform(Phit,pvec_)-1); 
[theta_con,MLval_con] = fmincon(MLfit,startp,A,b,A,b,0,100,@(pvec_)q_add_eq(pvec_,q_con));
if MLval<MLval_con && all(qform(Phit,params1)<=1)
    params = params1;
else
    params = theta_con;
end

% choose trials to predict
trials = trial_ix(groupProtocol,shelterSpace,groupSize,missing);
trEnd = zeros(1,size(Q1,1)); % we will store the time at which the trials ended
for j = 1:size(Q1,1) % iterate through each trial
    trEnd(j) = find(Q1(j,:)==Q1(j,end),1,'first');
end
nevac_b = zeros(numel(trials),1);
nevac_t = zeros(numel(trials),n_ts);
for tr = trials(:)'
    % randomly generate N individual evac times, p simulations
    % for comments see similar section above
    q = @(Phit_) qform(Phit_,params);
    tbins = 0:1:60;
    a = @(t_,P_,theta_) a_vec(t_,P_,theta_{1},theta_{2}); % this returns a prop. vector of only lower-triangular vals
    aix = tril(ones(N+1),-1);  % indices of lower-triangular entries
    theta = {makeA(N,q(Phit)); Q1(tr,:)}; % these are the non-explicitly-time-dependent input parameters
    P0 = [1;zeros(N,1)];   % this is the initial condition
    t0 = 0;                % starting time of simulation
    endx = find(Q1(tr,:)==Q1(tr,end),1,'first');
    T = endx + 1 - (endx==length(Q1(tr,:)));  % duration of simulation
    [ix,jx] = ind2sub(size(aix),find(~~aix));
    sx = sub2ind(size(aix),find(~~aix));
    nu = sparse(jx,1:length(sx),-1,N+1,length(sx)) + sparse(ix,1:length(sx),1);
    
    % this generates the samples with Gillespie algorithm...
    [t_mat,n_mat,r_mat,ns_mat,br_mat] = gillespie_tdep(theta,a,nu,P0,t0,T,tbins,p);
    
    % this creates a list of evac times & group IDs
    % randomly generate groups of size groupSize
    grpID1 = repmat(1:floor(N/groupSize),1,groupSize);
    grpID1 = [grpID1 zeros(1,N-length(grpID1))];
    grpIDs = zeros(p,N);
    n_evac = zeros(size(r_mat)); 
    n_evac(~~r_mat) = ix(r_mat(~~r_mat))-jx(r_mat(~~r_mat));
    t_evac = zeros(max(sum(n_evac)),p).*NaN;
    for px=1:p
      tx = [];
      for ne = removeval(unique(n_evac),0)'
         tx = [tx; repmat(t_mat(n_evac(:,px)==ne,px),ne,1)];
      end
      t_evac(1:sum(n_evac(:,px)),px) = sort(tx)';
      grpIDs(px,:) = grpID1(randperm(N));
    end
    
    % this creates cumulative evac plots given group IDs and evac times
    [Cgrp,tgrp,Cbin] = cum_evac(t_evac,grpIDs,groupProtocol,tbins);
    % plot cumulative evacuated
    figure; plot(C1(tr,:),'k--'); hold on; % observed in experiment
            plot(Q1(tr,:).*N,':');         % Phit trajectory
            for i=1:10%p
                plot(tgrp(:,i),Cgrp(:,i),'-o'); hold all;
                %plot(t_mat(r_mat(:,i)~=0,i),ix(r_mat(r_mat(:,i)~=0,i))-1,':o'); hold all;
            end
            plot(tbins,mean(Cbin,2),'-k');
            title(['trial ' num2str(tr) ', ' num2str(p) ' samples']); axis([0 60 0 50])
            legend('empirical data','Phit trajectory (scaled by N)',...
                   groupProtocol,...%'individual',
                   'location','northwest');
            xlabel('time'); ylabel('cumulative no. evacuated');
    
    % for each relevant trial, compute nevac_b and nevac(t)
    nevac_t(trials==tr,1:trEnd(tr)) = mean(Cbin(2:trEnd(tr)+1,:),2)' - C1(tr,1:trEnd(tr));
    nevac_b(trials==tr) = mean(Cbin(trEnd(tr)+1,:)) - C1(tr,end);            
            
end  % end loop over trials

% plot nevac_b for each game
horm = Q1(:,end);
hits = NaN*zeros(numel(trials),1); misses = hits;
hits(~~horm(trials)) = nevac_b(~~horm(trials)); 
misses(~horm(trials)) = nevac_b(~horm(trials));
figure; plot(nevac_b); hold on; 
        plot(zeros(numel(trials),1),'--k');
        plot(hits,'or'); plot(misses,'ok');
        xlabel('trial number'); set(gca,'XTickLabel',trials);
        ylabel('error in total number evacuated');
        title(['number evacuated prediction error, ' groupProtocol ', ss=' ...
               num2str(shelterSpace) ', grp=' num2str(groupSize)]);
           
% save from several groupProtocols with...
nevac_t_fTG_50_5 = nevac_t; 
nevac_b_fTG_50_5 = nevac_b;
nevac_t_lTG_50_5 = nevac_t;
nevac_b_lTG_50_5 = nevac_b;
nevac_t_mR_50_5 = nevac_t;
nevac_b_mR_50_5 = nevac_b;

% plot & compare
figure;
for gP = 3
    switch gP
        case 1, toPlot = nevac_t_fTG_50_5;
                mean_abs_b(1,1) = mean(abs(nevac_b_fTG_50_5));
                mean_abs_b(2,1) = mean(mean(abs(nevac_t_fTG_50_5)));
                plot(1:60,toPlot',':b'); hold on;
                plot(1:60,mean(toPlot),'-b');
        case 2, toPlot = nevac_t_lTG_50_5;
                mean_abs_b(1,2) = mean(abs(nevac_b_lTG_50_5));
                mean_abs_b(2,2) = mean(mean(abs(nevac_t_lTG_50_5)));
                plot(1:60,toPlot',':m'); hold on;
                plot(1:60,mean(toPlot),'-m');
        case 3, toPlot = nevac_t_mR_50_5;
                mean_abs_b(1,3) = mean(abs(nevac_b_mR_50_5));
                mean_abs_b(2,3) = mean(mean(abs(nevac_t_mR_50_5)));
                plot(1:60,toPlot'); hold on;
                plot(1:60,toPlot'); hold on;
                plot(1:60,toPlot,'-');
    end
end
title(['number evacuated prediction error, ss=' ...
               num2str(shelterSpace) ', grp=' num2str(groupSize)]);
xlabel('time'); ylabel('number evacuated error');
%legend('fTG','lTG','mR');

figure; bar(mean_abs_b');
        xlabel('group protocol'); set(gca,'XTickLabel',{'fTG';'lTG';'mR'});
        ylabel('mean absolute error');
        title(['total evacuations, mean absolute prediction error, ss=' ...
               num2str(shelterSpace) ', grp=' num2str(groupSize)]);
        legend('trial end error','mean over timesteps');
%% Making an analytical model for group protocol master equation
% we must translate individual probabilities to group ones.
% at any given time t, we have the individual-state prob. P[n,t] that n
% individuals will have evacuated at time t. we translate to the group
% probability P[ng,t] that exactly ng groups will have evacuated at t:
% 
% P[ng,t] = \sum_n (P[n,t]*Ptrans[groupProtocol])
%
% Ptrans(groupProtocol) is time-independent but takes into account the
% combinatorical probability for ng groups to have evacuated given n, and 
% group assignments selected at random from N individuals. it will also
% depend on the group size Gs, which itself is constant within a trial.
%
% the easiest is Ptrans[firstToGo] or Ptrans[fTG], which is the probability
% that n random selections without replacement from N individuals will be
% distributed in ng distinct groups (where groups are of preselected size
% Gs). each group with an evacuated individual in it will have evacuated.
% as such, Ptrans[fTG] is equal to P_group(n,ng,N,Gs).
tic;
Ptrans = grp_trans(N,Gs,groupProtocol);
toc;
qform = @(Phit_,pv_) pv_(1).*(Phit_.^pv_(3))./(pv_(2).^pv_(3)+Phit_.^pv_(3));
%qform = @(Phit_,pv_) pv_(1).*Phit_.^pv_(2);
Phit = 0:0.1:1;
q = @(Phit_) qhill(Phit_,L,k,n);
%[Ti, Pi] = solve_master_binom(N,q,Q1(tr,:),Phit);
Pgrp = Pi*Ptrans';



% next we try to 

%%  Bayesian learner that takes into acccount..........
%   the shelter space and others' decisions.
%   we need to formulate the likelihood function of getting the highest
%   score, GIVEN shelter space, Phit trajectory, and decision model.
%
%   if EVERYONE has model q(L,k,n) and ss=25... the probability of
%   successfully evacuating at timestep t = probability of personal decision
%   to evacuate * probability of shelter space being open = probability of
%   personal decision to evacuate * 
