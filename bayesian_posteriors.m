% this script will read in likelihoods from a particular bayesian run
% and compound them into posterior probabilities over time.

trials = 1:1:160;
ntrials = length(trials);
d = 2;
ss_list = [5,25,50];
prange = cell(d,1);
prange{1} = linspace(0.01,1,50);
prange{2} = linspace(3,10,36);
psize = cellfun(@numel,prange)';
tag = '_pl_ind_fine_err-20-13';

% set up saving struct
A = struct('likelihoods',cell(numel(ss_list),1),...
           'posteriors',cell(numel(ss_list),1),...
           'param_space',cell(numel(ss_list),1),...
           'norm_max_posterior',cell(numel(ss_list),1),...
           'n_opt_strategies',cell(numel(ss_list),1),...
           'opt_params',cell(numel(ss_list),1),...
           'opt_strategies',cell(numel(ss_list),1)...
           );
       
for ss=ss_list     
    likelihood = zeros([ntrials,psize]);    % likelihoods at each trial
    posteriors = zeros([ntrials+1,psize]);  % posterior probs at each trial
    maxprob = zeros(numel(trials)+1,1);     % maximum normalized posterior prob
    n_max_strategy = zeros(numel(trials)+1,1); % number of optimal strategies
    prob(1,:,:,:) = ones(psize);            % set flat prior (unnormalized)
    n_con = 0;                              % track number of constrained pts
    for px = 1:prod(psize)                  % choose one param combination
        pix = ind2sub_var_dim(psize,px,1);     % give the dimensional indices
        par = index_each_cell(prange,pix);     % give corresponding param val
        if any(q_con(par)>0)     % check this proposed strategy with constraint
          prob(1,pix{:}) = 0;
          n_con = n_con+1;
        end
    end
    maxprob(1) = 1/(prod(psize)-n_con);     % set initial normalized max prob.
    n_max_strategy(1) = prod(psize)-n_con;  % init: all valid strategies optimal
  % now, update for each trial
    for tr=trials
        disp(['ss = ' num2str(ss) ', tr = ' num2str(tr)]);
        load('likelihoods/LH_ss5_tr93_pl_ind_fine_err-20-13.mat'); % 'likelihoods'
        % update posteriors and save everything
        likelihood(trials==tr,:,:,:) = shiftdim(likelihoods,-1);
        prob(find(trials==tr)+1,:,:,:) = prob(trials==tr,:,:,:).*shiftdim(likelihoods,-1);
        norm = sum(sum(sum(prob(find(trials==tr)+1,:,:,:))));
        prob(find(trials==tr)+1,:,:,:) = prob(find(trials==tr)+1,:,:,:)./norm;
        maxprob(find(trials==tr)+1) = max(unique(prob(find(trials==tr)+1,:,:,:)))./norm;
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
    A(ss_list==ss).likelihoods = likelihood;
    A(ss_list==ss).posteriors = prob;
    A(ss_list==ss).param_space = prange;
    A(ss_list==ss).norm_max_posterior = maxprob;
    A(ss_list==ss).n_opt_strategies = n_max_strategy;
    A(ss_list==ss).opt_params = opt_params;
    A(ss_list==ss).opt_strategies = strategies;
    
end  % end loop over ss

%save(['bayesianLH' tag '.mat'],'A');