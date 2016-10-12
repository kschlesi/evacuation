%%%% Bayesian Learner of evacuation game, for arbitrary decision strategy
%%%% form (powerlaw, hill function, etc.)

%% First, load data.
 
bins = 0:0.1:1.1;  % this gives the bin EDGES (rounding DOWN)
trials = 1:160;
[H,J,Theta,AvgCumEvac,missing,Q1,T1,P1,C1] = load_evac_data(0,trials,bins);
[ntrials,n_ts] = size(Q1);
[N,~] = size(P1);

%% Next, create simple Bayesian learner.
% this updates in the order that trials were presented. use ALL trials.
trials = 1:1:ntrials; trials = trials(:)';

% loss matrix used to evaluate strategies
lossmat = loss_matrix(6,10,2,0);
                         
% decision model, fcn of Phit and parameters to fit
params = [1.0000;6.6423];                   % best fit params for powerlaw %%%%%%FIX
qform = @(Phit_,pv_) pv_(1).*Phit_.^pv_(2); % functional form of powerlaw
% qform = @(Phit_,pv_) pv_(1).*(Phit_.^pv_(3))./(pv_(2).^pv_(3)+Phit_.^pv_(3));
Phit = 0:0.1:1; % list of Phits at which q will be evaluated
q_con = @(pvec_)(qform(Phit,pvec_)-1); % constraint: q_con<=0

% set search area of parameter space
prange = cell(numel(params),1);
prange{1} = linspace(0.01,1,50);
prange{2} = linspace(3,10,36);
% prange{1} = linspace(0.6,2,5); 
% prange{2} = linspace(0.6,1.1,5);
% prange{3} = linspace(5,35,4);
psize = cellfun(@numel,prange)';

ss_list = N; ss = N;
A = struct('posteriors',cell(numel(ss_list),1),...
           'param_space',cell(numel(ss_list),1),...
           'norm_max_posterior',cell(numel(ss_list),1),...
           'n_opt_strategies',cell(numel(ss_list),1),...
           'opt_params',cell(numel(ss_list),1),...
           'opt_strategies',cell(numel(ss_list),1),...
           'likelihoods',cell(numel(ss_list),1)...
           );

% initialize matrices (normalizing 'prob' happens later)
prob = zeros([numel(trials)+1,psize]);  % posterior probs at each timestep
maxprob = zeros(numel(trials)+1,1);     % maximum normalized posterior prob
n_max_strategy = zeros(numel(trials)+1,1); % number of optimal strategies
prob(1,:,:,:) = ones(psize);            % set flat prior (unnormalized)
likes = zeros([numel(trials),psize]);   % raw likelihoods at each timestep
n_con = 0;                              % track number of constrained pts
for px = 1:prod(psize)                     % choose one param combination
    pix = ind2sub_var_dim(psize,px,1);     % give the dimensional indices
    par = index_each_cell(prange,pix);     % give corresponding param val
    if any(q_con(par)>0)     % check this proposed strategy with constraint
      prob(1,pix{:}) = 0;
      n_con = n_con+1;
    end
end
maxprob(1) = 1/(prod(psize)-n_con);     % set initial normalized max prob.
n_max_strategy(1) = prod(psize)-n_con;  % init: all valid strategies optimal

% iterate through each trial and update posterior prob of max score
for tr = trials
disp(['trial ' num2str(tr) ', ' num2str(find(trials==tr)) ' of ' num2str(numel(trials))]);
didHit = Q1(tr,end);
tic;
L_maxscore = zeros(psize);           % likelihood of achieving max score
for px=1:prod(psize)
    pix = ind2sub_var_dim(psize,px,1);     % give the dimensional indices
    par = index_each_cell(prange,pix);     % give corresponding param val
    % calculate evac likeliood given ss, Phit_traj, & params          
    [~,pEvac] = exp_score(qform,par,Phit,Q1(tr,:),lossmat,...
                          'AbsTol',1e-7,'RelTol',1e-8);
          if didHit % likelihood = pEvac
              L_maxscore(pix{:}) = pEvac;
          else      % likelihood = (1-pEvac)
              L_maxscore(pix{:}) = 1-pEvac;
          end
end
toc; % about 3 sec for 200 parameter space points (one trial)

% update posteriors, normalized max posterior, and # optimal strategies
likes(trials==tr,:,:,:) = shiftdim(L_maxscore,-1);
prob(find(trials==tr)+1,:,:,:) = prob(trials==tr,:,:,:).*shiftdim(L_maxscore,-1);
maxprob(find(trials==tr)+1) = max(unique(prob(find(trials==tr)+1,:,:,:)))./...
                                sum(sum(sum(prob(find(trials==tr)+1,:,:,:))));
n_max_strategy(find(trials==tr)+1) = numel(find(prob(find(trials==tr)+1,:,:,:)==max(unique(prob(find(trials==tr)+1,:,:,:)))));
end

% find strategy evolution
opt_params = zeros(sum(n_max_strategy==1),numel(psize)+1);
maxtrials = find(n_max_strategy==1);
for trx = 1:sum(n_max_strategy==1)
  tr = maxtrials(trx);  
  maxL = max(unique(prob(tr,:,:,:)));
  pix = ind2sub_var_dim(size(prob(tr,:,:,:)),find(prob(tr,:,:,:)==maxL),1);
  opt_params(trx,:) = [index_each_cell(prange,pix(2:end)),tr];
end
assert(opt_params(end,end)==numel(trials)+1);

% plot colorplot of strategy evolution
strategies = zeros(numel(maxtrials),numel(Phit));
for trx = 1:numel(maxtrials)
    tr = maxtrials(trx);
    strategies(trx,:) = qform(Phit,opt_params(trx,1:end-1));
end

% save Bayesian result
A.posteriors = prob;
A.param_space = prange;
A.norm_max_posterior = maxprob;
A.n_opt_strategies = n_max_strategy;
A.opt_params = opt_params;
A.opt_strategies = strategies;
A.likelihoods = likes;

%save('bayesian_pl_ss50_ind_fine_err-7-8.mat','A');

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
qform = @(Phit_,pv_) pv_(1).*Phit_.^pv_(2);
% qform = @(Phit_,pv_) pv_(1).*(Phit_.^pv_(3))./(pv_(2).^pv_(3)+Phit_.^pv_(3));
Phit = 0:0.1:1; % list of Phits at which q will be evaluated
q_con = @(pvec_)(qform(Phit,pvec_)-1); % constraint: q_con<=0

% set search area of parameter space
prange = cell(numel(params),1);
prange{1} = linspace(0.01,1,50);
prange{2} = linspace(3,10,36);
% prange{1} = linspace(0.6,2,5); 
% prange{2} = linspace(0.6,1.1,5);
% prange{3} = linspace(5,35,4);
psize = cellfun(@numel,prange)';

for ss=ss_list

% initialize matrices (normalizing 'prob' happens later)
prob = zeros([numel(trials)+1,psize]);  % posterior probs at each timestep
maxprob = zeros(numel(trials)+1,1);     % maximum normalized posterior prob
n_max_strategy = zeros(numel(trials)+1,1); % number of optimal strategies
prob(1,:,:,:) = ones(psize);            % set flat prior (unnormalized)
n_con = 0;                              % track number of constrained pts
A_bloc = zeros([N+1,N+1,numel(Phit),psize]);
for px = 1:prod(psize)                     % choose one param combination
    pix = ind2sub_var_dim(psize,px,1);     % give the dimensional indices
    par = index_each_cell(prange,pix);     % give corresponding param val
    if any(q_con(par)>0)     % check this proposed strategy with constraint
      prob(1,pix{:}) = 0;
      n_con = n_con+1;
    end
    Apar = makeA(N,qform(Phit,par),ss); % create A-mat for this strategy & ss
%     if ss<N     % update A-mat to disallow too many evacuations to shelter
%       for i=1:length(Phit)
%         zero blocked transitions and transfer probabilities to allowed ones
%         Apar(ss+1,:,i) = Apar(ss+1,:,i) + sum(Apar(ss+2:end,:,i));
%         Apar(ss+2:end,:,i) = 0; 
%         Apar(:,ss+2:end,i) = 0;
%         zero and re-compute diagonal
%         Apar(:,:,i) = Apar(:,:,i).*(~eye(size(Apar(:,:,i)))); assert(all(~diag(Apar(:,:,i))));
%         Apar(:,:,i) = Apar(:,:,i) - diag(sum(Apar(:,:,i))); %assert(all(~sum(A)));
%       end
%     end
    A_bloc(:,:,:,pix{:}) = Apar;
end
maxprob(1) = 1/(prod(psize)-n_con);     % set initial normalized max prob.
n_max_strategy(1) = prod(psize)-n_con;  % init: all valid strategies optimal

% iterate through each trial and update posterior prob of max score
for tr = trials
  disp(['trial ' num2str(tr) ', ' num2str(find(trials==tr)) ' of ' num2str(numel(trials))]);
  didHit = Q1(tr,end);
  %[gP,~,Gs] = trial_ix(tr);              % set group protocol, group size
  tic;
  L_maxscore = zeros(psize);              % likelihood of achieving max score
  for px=1:prod(psize)
    pix = ind2sub_var_dim(psize,px,1);     % give the dimensional indices
    par = index_each_cell(prange,pix);     % give corresponding param val
    % calculate evac likeliood given ss, Phit_traj, & params
    pEvac = succ_evac_prob(qform,par,Phit,tr,Q1(tr,:),lossmat,N,...
                           A_bloc(:,:,:,pix{:}),ss);
    if didHit % likelihood = pEvac
        L_maxscore(pix{:}) = pEvac;
    else      % likelihood = (1-pEvac)
        L_maxscore(pix{:}) = 1-pEvac;
    end
  end
  toc; % about 3 sec for 200 parameter space points (one trial)
  % update posteriors, normalized max posterior, and # optimal strategies
  prob(find(trials==tr)+1,:,:,:) = prob(trials==tr,:,:,:).*shiftdim(L_maxscore,-1);
  maxprob(find(trials==tr)+1) = max(unique(prob(find(trials==tr)+1,:,:,:)))./...
                                sum(sum(sum(prob(find(trials==tr)+1,:,:,:))));
  n_max_strategy(find(trials==tr)+1) = numel(find(prob(find(trials==tr)+1,:,:,:)==max(unique(prob(find(trials==tr)+1,:,:,:)))));
end

% find strategy evolution
opt_params = zeros(sum(n_max_strategy==1),numel(psize)+1);
maxtrials = find(n_max_strategy==1);
for trx = 1:sum(n_max_strategy==1)
  tr = maxtrials(trx);  
  maxL = max(unique(prob(tr,:,:,:)));
  pix = ind2sub_var_dim(size(prob(tr,:,:,:)),find(prob(tr,:,:,:)==maxL),1);
  opt_params(trx,:) = [index_each_cell(prange,pix(2:end)),tr];
end
assert(opt_params(end,end)==numel(trials)+1);

% plot colorplot of strategy evolution
strategies = zeros(numel(maxtrials),numel(Phit));
for trx = 1:numel(maxtrials)
    tr = maxtrials(trx);
    strategies(trx,:) = qform(Phit,opt_params(trx,1:end-1));
end

% save relevant things for this ss
A(ss_list==ss).posteriors = prob;
A(ss_list==ss).param_space = prange;
A(ss_list==ss).norm_max_posterior = maxprob;
A(ss_list==ss).n_opt_strategies = n_max_strategy;
A(ss_list==ss).opt_params = opt_params;
A(ss_list==ss).opt_strategies = strategies;
  
end  % end loop over ss

%save('bayesian_pl_ss50_ind_fine_err-40-13.mat','A');

%% what would the Bayesian player score on entire game?
lossmat = @(didHit_,didEvac_) 10*(didHit_.*~didEvac_) + ...
                             6*(didHit_.*didEvac_) + ...
                             0*(~didHit_.*~didEvac_) + ...
                             2*(~didHit_.*didEvac_) ;
qform = @(Phit_,pv_) pv_(1).*Phit_.^pv_(2);
%qform = @(Phit_,L_,k_,n_) L_.*(Phit_.^n_)./(k_.^n_+Phit_.^n_);
Phit = 0:0.1:1; 
q_con = @(pvec_)(qform(Phit,pvec_)-1); % constraint: q_con<=0
all_ss = zeros(numel(trials),1); ss_ix = zeros(numel(trials),1);
for tr=trials
    [~,all_ss(tr)] = trial_ix(tr);       % set shelter space capacity
end
ss_ix(all_ss==5) = 1; ss_ix(all_ss==25) = 2; ss_ix(all_ss==50) = 3;
load('bayesian_pl_ss_ind_fine_err-40-13.mat');

%%%%individual Bayesian player:
all_scr = zeros(numel(trials),1);
all_pEvac = zeros(numel(trials),1);
% choose random strategy to start
pxtry = num2cell(cellfun(@randi,num2cell(psize)));
ptry = index_each_cell(prange,pxtry);
while any(q_con(ptry)>0) % reject if strategy gives prob > 1
    pxtry = num2cell(cellfun(@randi,num2cell(psize)));
    ptry = index_each_cell(prange,pxtry);
end
for tr = trials
    % play trial tr
    didHit = Q1(tr,end);
    posts = A(ss_ix(tr)).posteriors(tr:tr+1,pxtry{:});
    pEvac = posts(2)/posts(1); if ~didHit; pEvac = 1-pEvac; end;
    all_scr(trials==tr) = -1*(lossmat(didHit,1)*pEvac + lossmat(didHit,0)*(1-pEvac));
    all_pEvac(trials==tr) = pEvac;
    if tr<trials(end)
      % update strategy for next trial tr+1
      ptry = A(ss_ix(tr+1)).opt_params(tr,1:end-1);
      pxtry = num2cell(find_each_cell(prange,ptry));
    end
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
lossmat = loss_matrix(6,10,2,0);
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

%% visualise landscapes
all_ss = zeros(numel(trials),1); ss_ix = zeros(numel(trials),1);
for tr=trials
    [~,all_ss(tr)] = trial_ix(tr);       % set shelter space capacity
end
ss_ix(all_ss==5) = 1; ss_ix(all_ss==25) = 2; ss_ix(all_ss==50) = 3;

% record movie
F(numel(trials)+1) = struct('cdata',[],'colormap',[]);
landscape = squeeze(A(1).posteriors(1,:,:,:));
figure;
surf(prange{2},prange{1},landscape);
F(1) = getframe(gcf);
for tr=trials
    %landscape = squeeze(A(ss_ix(tr)).posteriors(tr+1,:,:,:));
    %landscape = squeeze(A(3).posteriors(tr+1,:,:,:));
    landscape = squeeze(A.posteriors(tr+1,:,:,:));
    surf(prange{2},prange{1},landscape);
    F(tr+1) = getframe(gcf);
end

% play back movie
fig = figure;
movie(fig,F,1);
           
%% what would Bayesian strategy look like on individual games?
load('bayesian_pl_ss_ind_fine.mat');
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
