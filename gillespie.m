function [t_mat,n_mat,r_mat,ns_mat] = gillespie(theta,a,nu,n0,t0,T,nS)
% this is a discrete stochastic simulator of chemical systems, based on the
% gillespie-doob algorithm. it will run nS simulations of nR reactions,
% involving a total of nC chemical species.
%
% inputs: theta : nR x 1 vector of reaction parameters
%             a : nR x 1 vector of propensity functions (as function handles)
%            nu : nC x nR matrix denoting the net change in count of each
%                                species at the occurrence of each reaction
%            n0 : nC x 1 vector containing initial counts for all species
%            t0 : initial time of simulations
%             T : desired time duration of simulations
%            nS : number of simulations to run (optional, default 1)
%
% outputs: t_mat : max(ns_mat) x nS matrix of reaction times for each
%                                                                simulation
%          n_mat : nC x max(ns_mat) x nS matrix of population counts for
%                               each chemical species at each reaction time
%          r_mat : max(ns_mat) x nS matrix of reaction IDs
%         ns_mat : nS x 1 vector of total steps/reactions in each simulation

% set default value for nS
if nargin<7
    nS = 1;
end

% Gillespie algorithm initializations
t_mat = [repmat(t0,[1,nS]); zeros(T*max(theta),nS)];           % matrix to store reaction times
n_mat = [repmat(n0,[1,1,nS]), zeros(size(n0,1),T*max(theta),nS)];  % matrix to store n(t)
r_mat = [zeros([1,nS]);  zeros(T*max(theta),nS)];           % matrix to store reaction types
ns_mat = zeros(nS,1);

% loop over simulations
for sim=1:nS
    ns = 1;           % keeps track of number of steps run (inc. init cond)

    while t_mat(ns,sim)<(t0+T) && sum(a(n_mat(:,ns,sim),theta))>0
          % ensure that both time and reaction propensity remain

        % set t and n for this step
            t = t_mat(ns,sim);
            n = n_mat(:,ns,sim);
            disp(['run ' num2str(sim) ' of ' num2str(nS) ...
                  ', t = ' num2str(t) ' of ' num2str(t0+T)]);

        % draw next reaction time
            a0 = sum(a(n,theta));
            tnext = exprnd(1/a0);

        % draw next reaction
            % draw a uniform random number in [0,a0)
            r = a0;
            while r==a0 % reject if r==a0
                r = rand*a0;
            end
            % determine which reaction this indicates by proportionality in a
            % uniform interval (prob. of rxn j = aj/a0)
            rxn = find(cumsum(a(n,theta))>=r,1,'first');

        % update system with state-change matrix, increment time & step count
            n_mat(:,ns+1,sim) = n + nu(:,rxn);
            t_mat(ns+1,sim) = t + tnext;
            r_mat(ns+1,sim) = rxn;
            ns = ns+1;        

    end  % end algorithm while loop
    
    ns_mat(sim) = ns;  % store total number of reactions in simulation
    
end  % end loop over simulations

% cut off trailing zeros from storage matrices
n_mat = n_mat(:,1:max(ns_mat),:);
t_mat = t_mat(1:max(ns_mat),:);
r_mat = r_mat(1:max(ns_mat),:);

end