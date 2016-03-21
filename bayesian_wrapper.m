function outmessage = bayesian_wrapper(ss,tr,prange,psize,A_bloc,...
    qform,Phit,qcon,tolvec,tag,myDir)

% load Q1 and set game parameters
load([myDir '/evacuate_data.mat']);
Q1 = gameinfo;
[ntrials,~] = size(Q1);
trials = 1:1:ntrials;
N = 50;
lossmat = @(didHit_,didEvac_) 10*(didHit_.*~didEvac_) + ...
                              6*(didHit_.*didEvac_) + ...                           
                              0*(~didHit_.*~didEvac_) + ...
                              2*(~didHit_.*didEvac_) ;
                          
% compute likelihoods for this tr and ss
disp(['trial ' num2str(tr) ', ' num2str(find(trials==tr)) ' of ' num2str(numel(trials))]);
likelihoods = bayesian_game(psize,prange,qform,qcon,Phit,Q1,lossmat,N,...
                                                      A_bloc,tr,ss,tolvec);

% save likelihoods
save([myDir '/likelihoods/LH_ss' num2str(ss) '_tr' num2str(tr) tag '.mat']...
                                                        ,'likelihoods');
outmessage = ['likelihood success, ss = ' num2str(ss) ', tr = ' num2str(tr)];
                                                    
end