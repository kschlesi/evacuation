function pEvac = succ_evac_prob(qhill,params,Phit,tr,Phit_traj,lossmat,N,As,ss)

if nargin<8
    passAs = 0;
else
    passAs = 1;
end

[gP,nS,~] = trial_ix(tr);
if nargin<9
    ss=nS;
end

if all(strcmp(gP,'individual')) && all(ss==N)
    [~,pEvac] = exp_score(qhill,params,Phit,Phit_traj,lossmat);
else
  if passAs  
    [T,P] = solve_master_binom_ss(N,...
            @(Phit_) qhill(Phit_,params(1),params(2),params(3)),...
            Phit_traj,Phit,ss,As);
  else
    [T,P] = solve_master_binom_ss(N,...
          @(Phit_) qhill(Phit_,params(1),params(2),params(3)),...
          Phit_traj,Phit,ss);
  end
    pEvac = sum(bsxfun(@times,P,0:1:N),2)./N;  % for ind. at each timestep
    didHit = Phit_traj(end);
    trial_end_ix = find(T>=(find(Phit_traj==didHit,1,'first')+1),1,'first')-1;
    if ~numel(trial_end_ix)
        trial_end_ix = length(T);
    end
    pEvac = pEvac(trial_end_ix);    
end

end