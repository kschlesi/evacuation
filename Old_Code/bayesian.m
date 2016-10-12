
%% NOTE: This code is old (fits only for hill function) and an updated 
%% version can be found in bayesian_new
%%

%%%% Bayesian Learner of evacuation game.

%% First, load data.
 
bins = -0.05:0.1:1.05;  % this gives the bin EDGES
trials = 1:160;
[H,J,Theta,AvgCumEvac,missing,Q1,T1,P1,C1] = load_evac_data(0,trials,bins);
[ntrials,nts] = size(Q1);
[N,~] = size(P1);

%% Next, create simple Bayesian learner.
% this updates in the order that trials were presented. use ALL trials.
trials = 1:1:ntrials; trials = trials(:)';

% loss matrix used to evaluate strategies
lossmat = @(didHit_,didEvac_) 10*(didHit_.*~didEvac_) + ...
                             6*(didHit_.*didEvac_) + ...
                             0*(~didHit_.*~didEvac_) + ...
                             2*(~didHit_.*didEvac_) ;
                         
% decision model, fcn of Phit and parameters to fit
qhill = @(Phit_,L_,k_,n_) L_.*(Phit_.^n_)./(k_.^n_+Phit_.^n_);
Phit = 0:0.1:1; % list of Phits at which q will be evaluated
q_con = @(pvec_)(qhill(Phit,pvec_(1),pvec_(2),pvec_(3))-1); % constraint: q_con<=0

% set search area of parameter space
Lrange = linspace(0.6,2,5); nL = numel(Lrange);
krange = linspace(0.6,1.1,5); nk = numel(krange);
nrange = linspace(5,35,8); nn = numel(nrange);

ss_list=1;
A = struct('posteriors',cell(numel(ss_list),1),...
           'param_space',cell(numel(ss_list),1),...
           'norm_max_posterior',cell(numel(ss_list),1),...
           'n_opt_strategies',cell(numel(ss_list),1),...
           'opt_params',cell(numel(ss_list),1),...
           'opt_strategies',cell(numel(ss_list),1)...
           );

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

% save Bayesian result
A.posteriors = prob;
A.param_space = {Lrange,krange,nrange};
A.norm_max_posterior = maxprob;
A.n_opt_strategies = n_max_strategy;
A.opt_params = [Lt,kt,nt,trt];
A.opt_strategies = strategies;
  
save('bayesian_ss50_ind.mat','A');

%% Now, allow this Bayesian learner to take into account shelter space.
% what this means is that she will evaluate her own probability of
% evacuating SUCCESSFULLY rather than evacuating at all; she will modify
% her probability of evacuating by others' choices. (AND she will assume 
% others are doing the same?)
% this updates in the order that trials were presented. use ALL trials.
trials = 1:1:ntrials; trials = trials(:)';
ss_list = [5,25,50];
A = struct('posteriors',cell(numel(ss_list),1),...
           'param_space',cell(numel(ss_list),1),...
           'norm_max_posterior',cell(numel(ss_list),1),...
           'n_opt_strategies',cell(numel(ss_list),1),...
           'opt_params',cell(numel(ss_list),1),...
           'opt_strategies',cell(numel(ss_list),1)...
           );

% loss matrix used to evaluate strategies
lossmat = @(didHit_,didEvac_) 10*(didHit_.*~didEvac_) + ...
                             6*(didHit_.*didEvac_) + ...
                             0*(~didHit_.*~didEvac_) + ...
                             2*(~didHit_.*didEvac_) ;
                         
% decision model, fcn of Phit and parameters to fit
qhill = @(Phit_,L_,k_,n_) L_.*(Phit_.^n_)./(k_.^n_+Phit_.^n_);
Phit = 0:0.1:1; % list of Phits at which q will be evaluated
q_con = @(pvec_)(qhill(Phit,pvec_(1),pvec_(2),pvec_(3))-1); % constraint: q_con<=0

% set search area of parameter space
Lrange = linspace(0.6,2,5); nL = numel(Lrange);
krange = linspace(0.6,1.1,5); nk = numel(krange);
nrange = linspace(5,35,8); nn = numel(nrange);

for ss=ss_list

% initialize matrices (normalizing 'prob' happens later)
prob = zeros(nL,nk,nn,numel(trials)+1); % posterior probs at each timestep
maxprob = zeros(numel(trials)+1,1);     % maximum normalized posterior prob
n_max_strategy = zeros(numel(trials)+1,1); % number of optimal strategies
prob(:,:,:,1) = ones(nL,nk,nn);         % set flat prior (unnormalized)
n_con = 0;                              % track number of constrained pts
A_bloc = zeros(N+1,N+1,numel(Phit),nL,nk,nn);
for L=Lrange
  for k=krange
    for n=nrange 
    if any(q_con([L,k,n])>0) % check each proposed strategy with constraint
      prob(Lrange==L,krange==k,nrange==n,1) = 0;
      n_con = n_con+1;
    end
    ALkn = makeA(N,qhill(Phit,L,k,n));
    for i=1:length(Phit)
        % zero blocked transitions and transfer probabilities to allowed ones
        ALkn(ss+1,:,i) = ALkn(ss+1,:,i) + sum(ALkn(ss+2:end,:,i));
        ALkn(ss+2:end,:,i) = 0; 
        ALkn(:,ss+2:end,i) = 0;
        % zero and re-compute diagonal
        ALkn(:,:,i) = ALkn(:,:,i).*(~eye(size(ALkn(:,:,i)))); assert(all(~diag(ALkn(:,:,i))));
        ALkn(:,:,i) = ALkn(:,:,i) - diag(sum(ALkn(:,:,i))); %assert(all(~sum(A)));
    end
    A_bloc(:,:,:,Lrange==L,krange==k,nrange==n) = ALkn;
    end
  end
end
maxprob(1) = 1/(nn*nk*nL-n_con);        % set initial normalized max prob.
n_max_strategy(1) = nn*nk*nL-n_con;     % init: all valid strategies optimal

% iterate through each trial and update posterior prob of max score
for tr = trials
disp(['trial ' num2str(tr) ', ' num2str(find(trials==tr)) ' of ' num2str(numel(trials))]);
didHit = Q1(tr,end);
%[gP,ss,Gs] = trial_ix(tr);              % set group protocol, shelter, group size
tic;
L_maxscore = zeros(nL,nk,nn);           % likelihood of achieving max score
for L=Lrange
  for k=krange
    for n=nrange
          % calculate evac likeliood given ss, Phit_traj, & params
          pEvac = succ_evac_prob(qhill,[L,k,n],Phit,tr,Q1(tr,:),lossmat,N,...
                                 A_bloc(:,:,:,Lrange==L,krange==k,nrange==n),ss);
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
[Lox,kox,nox] = ind2sub(size(prob(:,:,:,end)),find(prob(:,:,:,end)==max(unique(prob(:,:,:,end)))));
Lopt = Lrange(Lox); kopt = krange(kox); nopt = nrange(nox);
Lt(end) = Lopt; kt(end) = kopt; nt(end) = nopt; 
trt(end) = trials(end);

% plot colorplot of strategy evolution
strategies = zeros(numel(trials),numel(Phit));
for tr = trials
    strategies(trials==tr,:) = qhill(Phit,Lt(trials==tr),kt(trials==tr),nt(trials==tr));
end

% save relevant things for this ss
A(ss_list==ss).posteriors = prob;
A(ss_list==ss).param_space = {Lrange,krange,nrange};
A(ss_list==ss).norm_max_posterior = maxprob;
A(ss_list==ss).n_opt_strategies = n_max_strategy;
A(ss_list==ss).opt_params = [Lt,kt,nt,trt];
A(ss_list==ss).opt_strategies = strategies;
  
end  % end loop over ss

%save('bayesian_pl_ss_ind.mat','A');

%% combine things....
for blah=blah1
    load('bayesian_ss5_ind.mat'); Ass5 = A;
    load('bayesian_ss25_ind.mat'); Ass25 = A;
    load('bayesian_ss50_ind.mat'); Ass50 = A;

    % check -- are 5s (and 50s) the same?
    load('bayesian_ssALL5_ind.mat');
    assert(all(A==Ass5));

    % make one large A and save
    A = struct('posteriors',cell(numel(ss_list),1),...
               'param_space',cell(numel(ss_list),1),...
               'norm_max_posterior',cell(numel(ss_list),1),...
               'n_opt_strategies',cell(numel(ss_list),1),...
               'opt_params',cell(numel(ss_list),1),...
               'opt_strategies',cell(numel(ss_list),1)...
               );
    A(1).posteriors = Ass5.posteriors;
    A(1).param_space = Ass5.param_space;
    A(1).norm_max_posterior = Ass5.norm_max_posterior;
    A(1).n_opt_strategies = Ass5.n_opt_strategies;
    A(1).opt_params = Ass5.opt_params;
    A(1).opt_strategies = Ass5.opt_strategies;

    A(2).posteriors = Ass25.posteriors;
    A(2).param_space = Ass25.param_space;
    A(2).norm_max_posterior = Ass25.norm_max_posterior;
    A(2).n_opt_strategies = Ass25.n_opt_strategies;
    A(2).opt_params = Ass25.opt_params;
    A(2).opt_strategies = Ass25.opt_strategies;

    A(3).posteriors = Ass50.posteriors;
    A(3).param_space = Ass50.param_space;
    A(3).norm_max_posterior = Ass50.norm_max_posterior;
    A(3).n_opt_strategies = Ass50.n_opt_strategies;
    A(3).opt_params = Ass50.opt_params;
    A(3).opt_strategies = Ass50.opt_strategies;

    %save('bayesian_ssALL_ind.mat','A');
end
%% what would the Bayesian player score on entire game?
lossmat = @(didHit_,didEvac_) 10*(didHit_.*~didEvac_) + ...
                             6*(didHit_.*didEvac_) + ...
                             0*(~didHit_.*~didEvac_) + ...
                             2*(~didHit_.*didEvac_) ;
qhill = @(Phit_,L_,k_,n_) L_.*(Phit_.^n_)./(k_.^n_+Phit_.^n_);
Phit = 0:0.1:1; 
q_con = @(pvec_)(qhill(Phit,pvec_(1),pvec_(2),pvec_(3))-1); % constraint: q_con<=0
all_ss = zeros(numel(trials),1); ss_ix = zeros(numel(trials),1);
for tr=trials
    [~,all_ss(tr)] = trial_ix(tr);       % set shelter space capacity
end
ss_ix(all_ss==5) = 1; ss_ix(all_ss==25) = 2; ss_ix(all_ss==50) = 3;
load('bayesian_ss_ind.mat');

%%%%individual Bayesian player:
all_scr = zeros(numel(trials),1);
all_pEvac = zeros(numel(trials),1);
% choose random strategy to start
ptry = [Lrange(randi(length(Lrange))),krange(randi(length(krange))),nrange(randi(length(nrange)))];
while any(q_con(ptry)>0) % reject if strategy gives prob > 1
    ptry = [Lrange(randi(length(Lrange))),krange(randi(length(krange))),nrange(randi(length(nrange)))];
end
for tr = trials
    % play trial
    didHit = Q1(tr,end);
    posts = A(ss_ix(tr)).posteriors(Lrange==ptry(1),krange==ptry(2),nrange==ptry(3),tr:tr+1);
    pEvac = posts(2)/posts(1); if ~didHit; pEvac = 1-pEvac; end;
    all_scr(trials==tr) = -1*(lossmat(didHit,1)*pEvac + lossmat(didHit,0)*(1-pEvac));
    all_pEvac(trials==tr) = pEvac;
    % update strategy
    ptry = [Lt(trials==tr),kt(trials==tr),nt(trials==tr)];
end

% bayesian expected losses per game
ix5 = zeros(size(all_scr))*NaN; ix25 = ix5; ix50 = ix5; 
hmask = ones(size(all_scr)); mmask = ones(size(all_scr));
hmask(Q1(:,end)==0) = NaN;
mmask(Q1(:,end)==1) = NaN;
ix5(all_ss==5) = all_scr(all_ss==5); 
ix25(all_ss==25) = all_scr(all_ss==25);
ix50(all_ss==50) = all_scr(all_ss==50);
figure; plot(ix5.*mmask,'or'); hold on; plot(ix5.*hmask,'*r');
        plot(ix25.*mmask,'ob'); plot(ix25.*hmask,'*b');
        plot(ix50.*mmask,'ok'); plot(ix50.*hmask,'*k');
        xlabel('trial'); ylabel('expected loss'); 
        title('Bayesian Performance (accounting for shelter space)');
        legend('miss, ss=5','hit, ss=5','miss, ss=25',...
               'hit, ss=25','miss, ss=50','hit, ss=50');
           
% bayesian prob. of evacuation per game
pe5 = zeros(size(all_scr))*NaN; pe25 = pe5; pe50 = pe5; 
pe5(all_ss==5) = all_pEvac(all_ss==5); 
pe25(all_ss==25) = all_pEvac(all_ss==25);
pe50(all_ss==50) = all_pEvac(all_ss==50);
figure; plot(pe5.*mmask,'or'); hold on; plot(pe5.*hmask,'*r');
        plot(pe25.*mmask,'ob'); plot(pe25.*hmask,'*b');
        plot(pe50.*mmask,'ok'); plot(pe50.*hmask,'*k');
        xlabel('trial'); ylabel('prob. of successful evacuation'); 
        title('Bayesian Performance (accounting for shelter space)');
        legend('miss, ss=5','hit, ss=5','miss, ss=25',...
               'hit, ss=25','miss, ss=50','hit, ss=50');

% figure; leg = cell(numel(find(all_scr<-9)),1);
% for tr = find(all_scr<-9)'
%     plot(Q1(tr,:)+0.0005*tr); hold on; 
%     leg{find(all_scr<-9)==tr} = ['trial ' num2str(tr)];
% end
% legend(leg,'location','northwest'); xlabel('time'); ylabel('Phit');
% title('No-Evacuate Disaster Hits for Bayesian Optimal Player');

%% comparison plots to actual behavior
lossmat = @(didHit_,didEvac_) 10*(didHit_.*~didEvac_) + ...
                               6*(didHit_.*didEvac_) + ...
                               0*(~didHit_.*~didEvac_) + ...
                               2*(~didHit_.*didEvac_) ;
didHit = Q1(:,end);
all_pEvac = C1(:,end)./N;
all_scr = -1*(lossmat(didHit,1).*all_pEvac + lossmat(didHit,0).*(1-all_pEvac));
all_pEvac(missing) = NaN; all_scr(missing) = NaN;

% actual average losses per game
ix5 = zeros(size(all_scr))*NaN; ix25 = ix5; ix50 = ix5; 
hmask = ones(size(all_scr)); mmask = ones(size(all_scr));
hmask(Q1(:,end)==0) = NaN;
mmask(Q1(:,end)==1) = NaN;
ix5(all_ss==5) = all_scr(all_ss==5); 
ix25(all_ss==25) = all_scr(all_ss==25);
ix50(all_ss==50) = all_scr(all_ss==50);
figure; plot(ix5.*mmask,'or'); hold on; plot(ix5.*hmask,'*r');
        plot(ix25.*mmask,'ob'); plot(ix25.*hmask,'*b');
        plot(ix50.*mmask,'ok'); plot(ix50.*hmask,'*k');
        xlabel('trial'); ylabel('average loss'); 
        title('Empirical Performance (accounting for shelter space)');
        legend('miss, ss=5','hit, ss=5','miss, ss=25',...
               'hit, ss=25','miss, ss=50','hit, ss=50');
           
% empirical evacuation prob. per game
pe5 = zeros(size(all_scr))*NaN; pe25 = pe5; pe50 = pe5; 
pe5(all_ss==5) = all_pEvac(all_ss==5); 
pe25(all_ss==25) = all_pEvac(all_ss==25);
pe50(all_ss==50) = all_pEvac(all_ss==50);
figure; plot(pe5.*mmask,'or'); hold on; plot(pe5.*hmask,'*r');
        plot(pe25.*mmask,'ob'); plot(pe25.*hmask,'*b');
        plot(pe50.*mmask,'ok'); plot(pe50.*hmask,'*k');
        xlabel('trial'); ylabel('empirical prob. of successful evacuation'); 
        title('Empirical Performance (accounting for shelter space)');
        legend('miss, ss=5','hit, ss=5','miss, ss=25',...
               'hit, ss=25','miss, ss=50','hit, ss=50');
           
%% what would Bayesian strategy look like on individual games?
load('bayesian_ss_ind.mat');
all_ss = zeros(numel(trials),1); ss_ix = zeros(numel(trials),1);
for tr=trials
    [~,all_ss(tr)] = trial_ix(tr);       % set shelter space capacity
end
ss_ix(all_ss==5) = 1; ss_ix(all_ss==25) = 2; ss_ix(all_ss==50) = 3;

% all ss=50 games
strategies = A(3).opt_strategies(all_ss==50,:);
figure; bcolor(strategies,Phit,trials(all_ss==50)); 
        c = colorbar; title('Bayesian strategy evolution, ss = 50');
        %set(gca,'XTick',1:5:ntrials,'XTickLabels',trials(all_ss==50));
        ylabel('trial'); xlabel('Phit value'); ylabel(c,'evacuation probability');
        
% all ss=25 games
strategies = A(2).opt_strategies(all_ss==25,:);
figure; bcolor(strategies,Phit,trials(all_ss==25)); 
        c = colorbar; title('Bayesian strategy evolution, ss = 25');
        %set(gca,'XTick',1:5:ntrials,'XTickLabels',trials(all_ss==25));
        ylabel('trial'); xlabel('Phit value'); ylabel(c,'evacuation probability');
        
% all ss=5 games
strategies = A(1).opt_strategies(all_ss==5,:);
figure; bcolor(strategies,Phit,trials(all_ss==5)); 
        c = colorbar; title('Bayesian strategy evolution, ss = 5');
        %set(gca,'XTick',1:5:ntrials,'XTickLabels',trials(all_ss==5));
        ylabel('trial'); xlabel('Phit value'); ylabel(c,'evacuation probability');
        
% all games in order
all_strategies = zeros(size(A(1).opt_strategies));
all_strategies(all_ss==5,:) = A(1).opt_strategies(all_ss==5,:);
all_strategies(all_ss==25,:) = A(2).opt_strategies(all_ss==25,:);
all_strategies(all_ss==50,:) = A(3).opt_strategies(all_ss==50,:);
figure; bcolor(all_strategies,Phit,trials); 
        c = colorbar; title('Bayesian strategy evolution, ss, all trials');
        %set(gca,'XTick',1:5:ntrials,'XTickLabels',trials(all_ss==5));
        ylabel('trial'); xlabel('Phit value'); ylabel(c,'evacuation probability');
        
% all unique strategies for each ss, colored by ss capacity
figure;
leg = cell(ntrials,1); lix = 0; nLeg = 0;
for ss=ss_list
    [u_params,trplot] = unique_rows(A(ss_list==ss).opt_params(:,1:3));
    nLeg = nLeg + numel(trplot);
    for tr=trplot'
        %leg{lix} = ['Final: \Lambda = ' num2str(Lopt) ', k = ' num2str(kopt) ', n = ' num2str(nopt)];
        lix = lix+1;
 fix       leg{lix} = ['Trial ' num2str(tr) ': \Lambda = ' num2str(Lt(tr)) ', k = ' num2str(kt(tr)) ', n = ' num2str(nt(tr))];
 fix       colorval = 1-find(all(bsxfun(@eq,unique_rows(strategies),qhill(Phit,Lt(tr),kt(tr),nt(tr)))'))/size(unique_rows(strategies),1);
            switch ss
                case 5, lcode = '--'; ccode = [0,colorval,colorval];
                case 25, lcode = ':'; ccode = [colorval,0,colorval];
                case 50, lcode = '.-'; ccode = [colorval,colorval,0];
            end
 fix       plot(Phit,all_strategies(tr,:),lcode,'Color',ccode);
    end
end    
xlabel('Phit'); ylabel('evacuation prob'); title('Bayesian strategy evolution');
legend(leg{1:nLeg},'location','northwest');
