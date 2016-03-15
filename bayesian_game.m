function L_maxscore = bayesian_game(psize,prange,qform,Phit,Q1,lossmat,N,...
                                    A_bloc,tr,ss)
% given all game info, this function computes the likelihood of successful
% evacuation for a single trial (tr) with a given shelter space (ss).
% required info: psize = dx1 vector, dimensions of strategy parameter space
%                prange = dx1 cell array, 
%                qform = function giving probability of individual evac
%                        based on Phit (nargin=2, arg1=Phit, arg2=dx1 param vector)
%                Phit = vector of possible Phit values
%                Q1 = (ntrials x time) matrix of disaster Phit trajectories
%                lossmat = function giving loss based on evacuation
%                          decision (nargin=2, arg1=didHit, arg2=didEvac)
%                N = scalar number of participants
%                A_bloc = possible generator matrices for markov master eqn
%                         for each Phit value and ss value
%                tr = scalar trial number
%                ss = scalar initial shelter space capacity
% output: L_maxscore, a matrix of size [psize] containing likelihood for
%                     successful evac at each parameter value in 'prange'

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

end