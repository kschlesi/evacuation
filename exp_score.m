function [scr,pEvac] = exp_score(q_model,pvec,Phit,Phit_traj,lossmat)
	% need nargin q_model = 2;
    % need nargin lossmat = 2;
	q = @(Phit_) q_model(Phit_,pvec);
	[T,P] = solve_master_binom(1,q,Phit_traj,Phit);
    didHit = Phit_traj(end);
    trial_end_ix = find(T>=(find(Phit_traj==didHit,1,'first')+1),1,'first')-1;
    if ~numel(trial_end_ix)
        trial_end_ix = length(T);
    end
    pEvac = P(trial_end_ix,end);
    scr = -1*(lossmat(didHit,1)*pEvac + lossmat(didHit,0)*(1-pEvac));
    if isnan(scr)
        disp([didHit,pEvac]);
    end
end
    
            