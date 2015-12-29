function scr = score(q_model,pvec,Phit,Phit_traj,lossmat)
	% need nargin q_model = 4;
    % need nargin lossmat = 2;
	q = @(Phit_) q_model(Phit_,pvec(1),pvec(2),pvec(3));
	[T,P] = solve_master_binom(1,q,Phit_traj,Phit);
    didHit = Phit_traj(end);
    trial_end_ix = find(T>=(find(Phit_traj==didHit,1,'first')+1),1,'first')-1;
    if ~numel(trial_end_ix)
        trial_end_ix = length(T);
    end
    didEvac = P(trial_end_ix,end)>=0.5;
    scr = -1*lossmat(didHit,didEvac);
end
    
            