function [T,P] = solve_master_binom(N,q,Phit_traj,Phit_range,varargin)
% assuming iid individual decisions, where q(Phit,params) = prob. of one
% individual to evacuate at a given Phit, this function uses the binomial
% distribution to construct the transition matrix and 'ode45' to solve 
% the master equation from a given initial evacuation state for a given
% trial realization (experimental Phit sequence).
% default Pinit = P0=1, all other states Pi have Pi=0
% default T_range = one solution at each new Phit value
%
% GENERAL MASTER EQUATION
% dP[n,t+dt] = SUM[i~=n:N](P[i,t]*P[n,t+dt|i,t] - P[n,t]*P[i,t+dt|n,t])
% If n < i then P[n, t + dt|i, t] = 0 "we cannot leave the shelter"
% we use P[n|i] = [N-i  n-i] q^(n-i) (1-q)^(N-n)

% Time Variable: define range of timesteps in which to give solutions
if any(ismemvar(varargin,'tRange'))
  T_range = varargin{find(ismemvar(varargin,'tRange'),1,'first')+1};
else
  T_range = [0,length(Phit_traj)-1];   % one solution at each solver timestep
end

% Initial probabilities of each state:
if any(ismemvar(varargin,'InitialState'))
  P_init = varargin{find(ismemvar(varargin,'InitialState'),1,'first')+1};
else
  P_init = zeros(N+1,1);
  P_init(1) = 1;
end

% A matrix at every value of Phit_range
As = makeA(N,q(Phit_range));

% Solver Options
if any(ismemvar(varargin,'AbsTol'))
    aTol = varargin{find(ismemvar(varargin,'AbsTol'),1,'first')+1};
    options = odeset('AbsTol',aTol);
end
if any(ismemvar(varargin,'RelTol'))
    rTol = varargin{find(ismemvar(varargin,'RelTol'),1,'first')+1};
    options = odeset('RelTol',rTol);
end

% Solver: solve the equation defined in 'odefunc' over the time range,
% starting with initial conditions P_init. this returns the variable T,
% a ts x 1 vector containing a list of times at which solutions are given,
% and P, a ts x length(P_init) vector containing the value of each entry of
% P(t) at each timestep in T.
[T, P] = ode45(@odefunc, T_range, P_init, options, Phit_traj, As);

% % Plot Results: create a plot of the trajectories at each value of P
% figure
% plot(T, P);  % plots all trajectories as a function of time
% title('Probability vs Time');
% ylabel('Probability of being in State n');
% xlabel('Time (seconds)');
% % create figure legend
% leg = cell(numel(P),1);
% for n=1:numel(P)
%     leg{n} = ['P' num2str(n-1) ' State' num2str(n-1)];
% end
% legend(leg,'Location','EastOutside');

end  % end main function

function dP = odefunc(t, P, Phit_traj, As)
% This function defines the differential equation that will be called in
% the main function. The first two arguments MUST be t (the timestep 
% calculated internally by ode45) and the variable being solved for (P).
% Other parameters (A) must be passed as arguments.
% The output should be the equivalent of dP/dt for one timestep.
A = As(:,:,floor(Phit_traj(floor(t)+1)*10+1));
dP = A*P;  
%disp(t);
end