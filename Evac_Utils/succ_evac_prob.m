function pEvac = succ_evac_prob(qhill,params,Phit,tr,Phit_traj,lossmat,N,varargin)

[gP,nS,~] = trial_ix(tr);
if any(ismemvar(varargin,'ss'))
    ss = varargin{find(ismemvar(varargin,'ss'),1,'first')+1};
else
    if any(ismemvar(varargin,'shelterSpace'))
        ss = varargin{find(ismemvar(varargin,'shelterSpace'),1,'first')+1};
    else
        ss = nS;
    end
end

tolArgs = cell(4,1);
if any(ismemvar(varargin,'relTol'))
    tolArgs{1} = 'relTol';
    tolArgs{2} = varargin{find(ismemvar(varargin,'relTol'),1,'first')+1};
end
if any(ismemvar(varargin,'absTol'))
    tolArgs{3} = 'absTol';
    tolArgs{4} = varargin{find(ismemvar(varargin,'absTol'),1,'first')+1};
end

if all(strcmp(gP,'individual')) && all(ss==N)
    [~,pEvac] = exp_score(qhill,params,Phit,Phit_traj,lossmat);
else
  if any(ismemvar(varargin,'As'))  
    As = varargin{find(ismemvar(varargin,'As'),1,'first')+1};
    [T,P] = solve_master_binom_ss(N,...
            @(Phit_) qhill(Phit_,params(1),params(2),params(3)),...
            Phit_traj,Phit,'ss',ss,'As',As,tolArgs{:});
  else
    [T,P] = solve_master_binom_ss(N,...
          @(Phit_) qhill(Phit_,params(1),params(2),params(3)),...
          Phit_traj,Phit,'ss',ss,tolArgs{:});
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