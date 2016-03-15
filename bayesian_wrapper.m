function outmessage = bayesian_wrapper(ss,tr,prange,psize,A_bloc,myDir,tag)

% load Q1 and set game parameters
load([myDir '\evacuate_data.mat']);
Q1 = gameinfo;
[ntrials,~] = size(Q1);
trials = 1:1:ntrials;
N = 50;
lossmat = @(didHit_,didEvac_) 10*(didHit_.*~didEvac_) + ...
                              6*(didHit_.*didEvac_) + ...
                              0*(~didHit_.*~didEvac_) + ...
                              2*(~didHit_.*didEvac_) ;
                          
% set bayesian strategy and parameter space
qform = @(Phit_,pv_) pv_(1).*Phit_.^pv_(2);
Phit = 0:0.1:1; % list of Phits at which q will be evaluated
q_con = @(pvec_)(qform(Phit,pvec_)-1); % constraint: q_con<=0

% compute likelihoods for this tr and ss
disp(['trial ' num2str(tr) ', ' num2str(find(trials==tr)) ' of ' num2str(numel(trials))]);
likelihoods = bayesian_game(psize,prange,qform,Phit,Q1,lossmat,N,A_bloc,ss);
outmessage = ['likelihood success, ss = ' num2str(ss) ', tr = ' num2str(tr)];

% save likelihoods
save([myDir '\likelihoods\LH' tag '.mat'],'likelihoods');

end