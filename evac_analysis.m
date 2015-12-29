%%% script for running things with evacuation data

% load_evac_data 
bins = -0.05:0.1:1.05;  % this gives the bin EDGES
trials = [19,28,29,49,64,67,76,102,108,112,113,123,125,130,143,144];
[H,J,Theta,AvgCumEvac,missing,Q1,T1,P1,C1] = load_evac_data(0,trials,bins);
[ntrials,nts] = size(Q1);
[N,~] = size(P1);

% functionalities:
% try leaving out 0
%H = H(:,2:end); J = J(:,2:end);
%% we can (a) fit experimental data to a particular q function 
% with maximum likelihood estimation
qhill = @(Phit_,L_,k_,n_) L_.*(Phit_.^n_)./(k_.^n_+Phit_.^n_);
Phit = 0:0.1:1; %Phit = Phit(2:end);
qfit = @(L_,k_,n_) qhill(Phit(:,2:end),L_,k_,n_);
MLfit = @(pvec) -1*sum(((sum(H(:,2:end))-sum(J(:,2:end))).*log(1-qfit(pvec(1),pvec(2),pvec(3))) + sum(J(:,2:end)).*log(qfit(pvec(1),pvec(2),pvec(3)))));
%MLfitC = @(pvec) -1*sum(((cumuH(2:end)'-cumuJ(2:end)').*log(1-qfit(pvec(1),pvec(2),pvec(3))) + cumuJ(2:end)'.*log(qfit(pvec(1),pvec(2),pvec(3)))));
startp = [1;0.5;10];
options = optimoptions(@fminunc,'MaxFunEvals',10000);
[params1,MLval] = fminunc(MLfit,startp,options);
A = zeros(length(startp));
b = [0 0 0];
% constrains all prob. as a function of Phit to be <=1
q_con = @(pvec_)(qhill(Phit,pvec_(1),pvec_(2),pvec_(3))-1); 
[theta_con,MLval_con] = fmincon(MLfit,startp,A,b,A,b,0,100,@(pvec_)q_add_eq(pvec_,q_con));


if MLval<MLval_con || all(qhill(Phit,params1(1),params(2),params(3))<=1)
    params = params1;
else
    params = theta_con;
end

figure; plot(Phit(2:end),qfit(params(1),params(2),params(3))); hold on;
legend(['\Lambda = ' num2str(params(1)) ', k = ' num2str(params(2)) ', n = ' num2str(params(3))]);
xlabel('Phit'); ylabel('evacuation probability'); title('best-fit hill model, ind. trials, shelterSpace=50');
        
%% we can (b) solve the master equation for a given trial and q function
%for tr=trials
tr=19;
L = params(1);
k = params(2);
n = params(3);
Phit = 0:0.1:1;
q = @(Phit_) qhill(Phit_,L,k,n);
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
L = params(1);
k = params(2);
n = params(3);
q = @(Phit_) qhill(Phit_,L,k,n);
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
L = params(1);
k = params(2);
n = params(3);
tbins = 0:1:60;
P0_ind = 0;            % initial value of "evacuated" state
t0 = 0; T = 60;        % starting time and duration of simulation
q_ind = @(Phit_) qhill(Phit_,L,k,n);
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
tr = 19;
L = params(1);
k = params(2);
n = params(3);
Phit = 0:0.1:1;
q = @(Phit_) qhill(Phit_,L,k,n);
[T, P] = solve_master_binom(1,q,Q1(tr,:),Phit);

figure; plot(T,P(:,end)'); hold on; 
        plot(1:1:60,Q1(tr,:)); 
        legend('prob. of evacuation','Phit','location','northwest');
        xlabel('time'); title(['Individual Evacuation Probability, Trial ' num2str(tr)]);
        
%% fit for STATIC OPTIMAL behavior
% for a given set of trials ('trials')
[Qu,trials] = unique_rows(Q1);
lossmat = @(didHit,didEvac) 10*(didHit.*~didEvac) + ...
                             6*(didHit.*didEvac) + ...
                             0*(~didHit.*~didEvac) + ...
                             2*(~didHit.*didEvac) ;
qhill = @(Phit_,L_,k_,n_) L_.*(Phit_.^n_)./(k_.^n_+Phit_.^n_);
Phit = 0:0.1:1; %Phit = Phit(2:end);
% for a single trial: as an example
ptry = params;
scr = exp_score(qhill,ptry,Phit,Q1(tr,:),lossmat);

% for all trials:
[tscr,scrs] = total_score(qhill,ptry,Phit,trials,Q1,lossmat);

% doing the fit:
%startp = params;
%startp = sopt_params_con;
startp = [Lrange(Lix),krange(kix),nrange(nix)];
SFit = @(ptry_) -1*total_score(qhill,ptry_,Phit,trials,Q1,lossmat);
%options = optimoptions(@fminunc,'MaxFunEvals',10000);
%[sopt_params,sopt_score] = fminunc(SFit,startp,options);
A = zeros(length(startp));
b = [0 0 0];
% constrains all prob. as a function of Phit to be <=1
q_con = @(pvec_)(qhill(Phit,pvec_(1),pvec_(2),pvec_(3))-1); 
[sopt_params_con,sopt_score_con] = fmincon(SFit,startp,A,b,A,b,0,100,@(pvec_)q_add_eq(pvec_,q_con));
%params = sopt_params;

% ok, mapping over the range:
Lrange = linspace(0.6,2,5);
krange = linspace(0.6,1.1,5);
nrange = linspace(5,35,8);
scoremat = zeros(numel(Lrange),numel(krange),numel(nrange));
tic;
for L=Lrange
  for k=krange
    for n=nrange
      if any(q_con([L,k,n])>0)
          scoremat(Lrange==L,krange==k,nrange==n) = NaN;
      else
          scoremat(Lrange==L,krange==k,nrange==n) = SFit([L,k,n]);
      end
    end
  end
end
toc; % about 47 sec for 200 points
best = min(min(min(scoremat)));
allscores = unique(scoremat);
[Lix,kix,nix] = ind2sub(size(scoremat),find(scoremat==best));
[Lix2,kix2,nix2] = ind2sub(size(scoremat),find((scoremat>allscores(1)).*(scoremat<allscores(15))));
[Lix3,kix3,nix3] = ind2sub(size(scoremat),find((scoremat>=allscores(15)).*(scoremat<allscores(35))));
[Lix4,kix4,nix4] = ind2sub(size(scoremat),find((scoremat>=allscores(35)).*(scoremat<allscores(55))));
[Lix5,kix5,nix5] = ind2sub(size(scoremat),find((scoremat>=allscores(55)).*(scoremat<allscores(75))));
% plot figures of all scores and max scores
figure; scatter3(Lrange(Lix),krange(kix),nrange(nix),'r'); hold on;
        scatter3(Lrange(Lix2),krange(kix2),nrange(nix2),'b');
        scatter3(Lrange(Lix3),krange(kix3),nrange(nix3),'k');
        scatter3(Lrange(Lix4),krange(kix4),nrange(nix4),'g');
        scatter3(Lrange(Lix5),krange(kix5),nrange(nix5),'y');
        title('best score parameters');
        xlabel('L'); ylabel('k'); zlabel('n');
        set(gca,'XTick',Lrange,'YTick',krange,'ZTick',nrange);
        axis([Lrange(1),Lrange(end),krange(1),krange(end),nrange(1),nrange(end)])
        legend([num2str(best) ', max score'],num2str(allscores(15)),...
                            num2str(allscores(35)),num2str(allscores(55)),...
                            num2str(allscores(75)),'location','northeast');
        
% what were the mistakes on optimal scores?
bx = [Lix,kix,nix];
all_opt = zeros(size(bx,1),numel(trials));
for i=1:size(bx,1)
    sumall = scoremat(bx(i,1),bx(i,2),bx(i,3));
    ptry = [Lrange(bx(i,1)),krange(bx(i,2)),nrange(bx(i,3))];
    [~,all] = total_score(qhill,ptry,Phit,trials,Q1,lossmat);
    all_opt(i,:) = all;
end
figure; bcolor([Q1(trials,end)';all_opt]); colorbar;

% plot maximum score found by mesh search of param space
figure;
%leg = cell(size(bx,1),1);
%for i=1:size(bx,1)
leg = cell(numel(find(krange(bx(:,2))~=0)),1);
for i = find(krange(bx(:,2))~=0)
    plot(Phit,qhill(Phit,Lrange(bx(i,1)),krange(bx(i,2)),nrange(bx(i,3)))); hold all; 
    leg{find(krange(bx(:,2))~=0)==i} = ['\Lambda = ' num2str(Lrange(bx(i,1))) ', k = ' num2str(krange(bx(i,2))) ', n = ' num2str(nrange(bx(i,3)))];
end
xlabel('Phit'); ylabel('Evacuation Probability'); legend(leg,'location','northwest');
title(['Static Optimal Evacuation Strategy (Hill Model), <score> = ' num2str(scoremat(Lix,kix,nix))]);

% plot maximum score found overall
figure;
    plot(Phit,qhill(Phit,sopt_params_con(1),sopt_params_con(2),sopt_params_con(3))); hold all; 
    leg = ['\Lambda = ' num2str(sopt_params_con(1)) ', k = ' num2str(sopt_params_con(2)) ', n = ' num2str(sopt_params_con(3))];
xlabel('Phit'); ylabel('Evacuation Probability'); legend(leg,'location','northwest');
title(['Static Optimal Evacuation Strategy (Hill Model), <score> = ' num2str(sopt_score_con)]);

%% fit for BAYESIAN OPTIMAL behavior
% this updates in the order that trials were presented. use ALL trials.
trials = 1:1:size(Q1,1); trials = trials(:)';
%trials = [19,28,29,49,64,67,76,102,108,112,113,123,125,130,143,144];

% loss matrix used to evaluate strategies
lossmat = @(didHit,didEvac) 10*(didHit.*~didEvac) + ...
                             6*(didHit.*didEvac) + ...
                             0*(~didHit.*~didEvac) + ...
                             2*(~didHit.*didEvac) ;
                         
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
tic;
L_maxscore = zeros(nL,nk,nn);           % likelihood of achieving max score
for L=Lrange
  for k=krange
    for n=nrange
          [~,pEvac] = exp_score(qhill,[L,k,n],Phit,Q1(tr,:),lossmat);  
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


%%%% to do still: make plots of highest probability areas over time

%% HOW does the Bayesian player score? 
lossmat = @(didHit,didEvac) 10*(didHit.*~didEvac) + ...
                             6*(didHit.*didEvac) + ...
                             0*(~didHit.*~didEvac) + ...
                             2*(~didHit.*didEvac) ;
qhill = @(Phit_,L_,k_,n_) L_.*(Phit_.^n_)./(k_.^n_+Phit_.^n_);
Phit = 0:0.1:1; 
q_con = @(pvec_)(qhill(Phit,pvec_(1),pvec_(2),pvec_(3))-1); % constraint: q_con<=0

%%%%individual Bayesian player:
all_scr = zeros(numel(trials),1);
% choose random strategy to start
ptry = [Lrange(randi(length(Lrange))),krange(randi(length(krange))),nrange(randi(length(nrange)))];
while any(q_con(ptry)>0) % reject if strategy gives prob > 1
    ptry = [Lrange(randi(length(Lrange))),krange(randi(length(krange))),nrange(randi(length(nrange)))];
end
for tr = trials
    % play trial
    all_scr(trials==tr) = exp_score(qhill,ptry,Phit,Q1(tr,:),lossmat);
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
q_sopt = @(Phit_) qhill(Phit_,sopt_params_con(1),sopt_params_con(2),sopt_params_con(3));

% choose random starting strategy from flat prior
ptry = [Lrange(randi(length(Lrange))),krange(randi(length(krange))),nrange(randi(length(nrange)))];
while any(q_con(ptry)>0) % reject if strategy gives prob > 1
    ptry = [Lrange(randi(length(Lrange))),krange(randi(length(krange))),nrange(randi(length(nrange)))];
end
for tr=trials'
    qb = @(Phit_) qhill(Phit_,ptry(1),ptry(2),ptry(3));
    [T, P] = solve_master_binom(N,qb,Q1(tr,:),Phit);
    Cexp = sum(bsxfun(@times,P,0:1:N),2);
    q = @(Phit_) qhill(Phit_,params(1),params(2),params(3));
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
bins = -0.05:0.1:1.05; [H,J] = load_evac_data(0,trials,bins);
qhill = @(Phit_,L_,k_,n_) L_.*(Phit_.^n_)./(k_.^n_+Phit_.^n_);
Phit = 0:0.1:1; qfit = @(L_,k_,n_) qhill(Phit(:,2:end),L_,k_,n_);
MLfit = @(pvec) -1*sum(((sum(H(:,2:end))-sum(J(:,2:end))).*log(1-qfit(pvec(1),pvec(2),pvec(3))) + sum(J(:,2:end)).*log(qfit(pvec(1),pvec(2),pvec(3)))));
startp = [1;0.5;10];
options = optimoptions(@fminunc,'MaxFunEvals',10000);
[params1,MLval] = fminunc(MLfit,startp,options);
A = zeros(length(startp)); b = [0 0 0]; q_con = @(pvec_)(qhill(Phit,pvec_(1),pvec_(2),pvec_(3))-1); 
[theta_con,MLval_con] = fmincon(MLfit,startp,A,b,A,b,0,100,@(pvec_)q_add_eq(pvec_,q_con));
if MLval<MLval_con || all(qhill(Phit,params1(1),params(2),params(3))<=1)
    params = params1;
else
    params = theta_con;
end

% choose trials to predict
trials = trial_ix(groupProtocol,shelterSpace,groupSize,missing);
for tr = trials(:)'
    % randomly generate N individual evac times, p simulations
    % for comments see similar section above
    q = @(Phit_) qhill(Phit_,params(1),params(2),params(3));
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
end

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
L = params(1);
k = params(2);
n = params(3);
qhill = @(Phit_,L_,k_,n_) L_.*(Phit_.^n_)./(k_.^n_+Phit_.^n_);
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