function [t_mat,n_mat,r_mat,ns_mat,br_mat] = gillespie_tdep(theta,a,nu,n0,t0,T,tbins,nS)
% this is a discrete stochastic simulator of chemical systems, based on the
% gillespie-doob algorithm. it will run nS simulations of nR reactions,
% involving a total of nC chemical species.
%
% inputs: theta : nR x 1 vector of reaction parameters
%             a : nR x 1 vector of propensity functions (as function handles)
%                *here, a can be time-dependent a(t,n,theta); specify tbins
%            nu : nC x nR matrix denoting the net change in count of each
%                                species at the occurrence of each reaction
%            n0 : nC x 1 vector containing initial counts for all species
%            t0 : initial time of simulations
%             T : desired time duration of simulations
%         tbins : times at which a changes (or edges of bins of constant a)
%            nS : number of simulations to run (optional, default 1)
%
% outputs: t_mat : max(ns_mat) x nS matrix of reaction times for each
%                                                                simulation
%          n_mat : nC x max(ns_mat) x nS matrix of population counts for
%                               each chemical species at each reaction time
%          r_mat : max(ns_mat) x nS matrix of reaction IDs
%         ns_mat : nS x 1 vector of total steps/reactions in each simulation

% set default value for nS
if nargin<8
    nS = 1;
end

% if a NOT time-dependent, call gillespie.m
if nargin(a)==2 || (all(~tbins) || ~numel(tbins))
    warning('using time-independent gillespie algorithm');
    br_mat = [];
    [t_mat,n_mat,r_mat,ns_mat] = gillespie(theta,a,nu,n0,t0,T,nS);
    return
end

% ensure tbins covers simulation time
tbins = tbins(:)';
try assert(tbins(1)<=t0)
catch 
    warning(['assuming constant propensity before t = ' num2str(tbins(1))]);
    tbins = [t0 tbins];
end
try assert(tbins(end)>=t0+T)
catch 
    warning(['assuming constant propensity after t = ' num2str(tbins(end))]);
    tbins = [tbins t0+T];
end

% define bin-finder function (rounds to floor of bin)
bin = @(t_) find([tbins,t_+tbins(end)]>t_,1,'first')-1;

% Gillespie algorithm initializations
t_mat = [repmat(t0,[1,nS]); zeros(T,nS)];           % matrix to store reaction times
n_mat = [repmat(n0,[1,1,nS]), zeros(size(n0,1),T,nS)];  % matrix to store n(t)
r_mat = [zeros([1,nS]); zeros(T,nS)];           % matrix to store reaction types
ns_mat = zeros(nS,1);
br_mat = zeros(nS,1);

% loop over simulations
for sim=1:nS
    ns = 1;           % keeps track of number of steps run (inc. init cond)

    while t_mat(ns,sim)<(t0+T) && sum(a(t_mat(ns,sim),n_mat(:,ns,sim),theta))>0
          % ensure that both time and reaction propensity remain

        % set t and n for this step
            t = t_mat(ns,sim);
            n = n_mat(:,ns,sim);
            disp(['run ' num2str(sim) ' of ' num2str(nS) ...
                  ', t = ' num2str(t) ' of ' num2str(t0+T)]);

        % draw next reaction time
            a0 = sum(a(t,n,theta));
            tnext = exprnd(1/a0);
            % ensure reaction occurs in timebin of constant propensity
            if bin(t)~=bin(t+tnext)
                %warning('timestep breaks discreteness');
                % advance to next time bin and update step count,
                % but leave state unchanged
                br_mat(sim) = br_mat(sim)+1;
                n_mat(:,ns+1,sim) = n;
                t_mat(ns+1,sim) = tbins(bin(t)+1);
                r_mat(ns+1,sim) = 0;
                ns = ns+1;
                continue
            end

        % draw next reaction
            % draw a uniform random number in [0,a0)
            r = a0;
            while r==a0 % reject if r==a0
                r = rand*a0;
            end
            % determine which reaction this indicates by proportionality in a
            % uniform interval (prob. of rxn j = aj/a0)
            rxn = find(cumsum(a(t,n,theta))>=r,1,'first');

        % update system with state-change matrix, increment time & step count
            n_mat(:,ns+1,sim) = n + nu(:,rxn);
            t_mat(ns+1,sim) = t + tnext;
            r_mat(ns+1,sim) = rxn;
            ns = ns+1;        

    end  % end algorithm while loop
    
    ns_mat(sim) = ns;  % store total number of reactions in simulation
    %disp(['condition 1: t = ' num2str(t_mat(ns,sim)) ', t_end = ' num2str(t0+T)]);
    %disp(['condition 2: sum of propensities = ' num2str(sum(a(t_mat(ns,sim),n_mat(:,ns,sim),theta)))]);
        
end  % end loop over simulations

% cut off trailing zeros from storage matrices
n_mat = n_mat(:,1:max(ns_mat),:);
t_mat = t_mat(1:max(ns_mat),:);
r_mat = r_mat(1:max(ns_mat),:);

end