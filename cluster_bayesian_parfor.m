% A script using the default parallel configuration to assign 
%
%  Using matlabpool and parfor
%
tic;
matlabpool OPEN 40;
myDir = pwd;

% load overall data, set game parameters
ntrials = 160;
N = 50;
ss_list = [5,25,50];
aT = 1e-20;
rT = 1e-13;
prange = cell(numel(params),1);
prange{1} = linspace(0.01,1,50);
prange{2} = linspace(3,10,36);
psize = cellfun(@numel,prange)';
tag = '_pl_ind_fine_err-20-13';

% create A_bloc for each ss
for ss=ss_list
    eval(['A_bloc' num2str(ss) ' = zeros([N+1,N+1,numel(Phit),psize]);']);
    for px=1:prod(psize)
    Apar = makeA(N,qform(Phit,par)); % create A-mat for this strategy & ss
    if ss<N     % update A-mat to disallow too many evacuations to shelter
      for i=1:length(Phit)
        % zero blocked transitions and transfer probabilities to allowed ones
        Apar(ss+1,:,i) = Apar(ss+1,:,i) + sum(Apar(ss+2:end,:,i));
        Apar(ss+2:end,:,i) = 0; 
        Apar(:,ss+2:end,i) = 0;
        % zero and re-compute diagonal
        Apar(:,:,i) = Apar(:,:,i).*(~eye(size(Apar(:,:,i)))); assert(all(~diag(Apar(:,:,i))));
        Apar(:,:,i) = Apar(:,:,i) - diag(sum(Apar(:,:,i))); %assert(all(~sum(A)));
      end
    end
    eval(['A_bloc' num2str(ss) '(:,:,:,pix{:}) = Apar;']);
    end 
end

% set number of tasks and pass to parallel loop
numtasks = ntrials*numel(ss_list);
outputarray = cell(numtasks,1);
parfor tasknum = 1:numtasks
    ss_list = [5,25,50];
    ssx = floor((tasknum-1)/ntrials)+1;
    ss = ss_list(ssx);
    switch ss
        case 5,  A_bloc = A_bloc5;
        case 25, A_bloc = A_bloc25;
        case 50, A_bloc = A_bloc50;
    end
    tr = mod(tasknum-1,ntrials)+1;
    outputarray{tasknum} = bayesian_wrapper(ss,tr,prange,psize,A_bloc,...
                                                        [aT,rT],tag,myDir);
end

% trials = 1:1:ntrials
% % set up saving struct
% A = struct('likelihoods',cell(numel(ss_list),1),...
%            'posteriors',cell(numel(ss_list),1),...
%            'param_space',cell(numel(ss_list),1),...
%            'norm_max_posterior',cell(numel(ss_list),1),...
%            'n_opt_strategies',cell(numel(ss_list),1),...
%            'opt_params',cell(numel(ss_list),1),...
%            'opt_strategies',cell(numel(ss_list),1)...
%            );
% posteriors = zeros([ntrials+1,psize]);  % posterior probs at each timestep
% maxprob = zeros(numel(trials)+1,1);     % maximum normalized posterior prob
% n_max_strategy = zeros(numel(trials)+1,1); % number of optimal strategies
% prob(1,:,:,:) = ones(psize);            % set flat prior (unnormalized)
% n_con = 0;                              % track number of constrained pts
% for px = 1:prod(psize)                  % choose one param combination
%     pix = ind2sub_var_dim(psize,px,1);     % give the dimensional indices
%     par = index_each_cell(prange,pix);     % give corresponding param val
%     if any(q_con(par)>0)     % check this proposed strategy with constraint
%       prob(1,pix{:}) = 0;
%       n_con = n_con+1;
%     end
% end
% maxprob(1) = 1/(prod(psize)-n_con);     % set initial normalized max prob.
% n_max_strategy(1) = prod(psize)-n_con;  % init: all valid strategies optimal
% for tr=trials
% % update posteriors and save everything
% prob(find(trials==tr)+1,:,:,:) = prob(trials==tr,:,:,:).*shiftdim(L_maxscore,-1);
% maxprob(find(trials==tr)+1) = max(unique(prob(find(trials==tr)+1,:,:,:)))./...
%                                 sum(sum(sum(prob(find(trials==tr)+1,:,:,:))));
% n_max_strategy(find(trials==tr)+1) = numel(find(prob(find(trials==tr)+1,:,:,:)==max(unique(prob(find(trials==tr)+1,:,:,:)))));
% end

matlabpool CLOSE;
for tasknum = 1:numtasks
    disp(outputarray{tasknum});
end
toc;
