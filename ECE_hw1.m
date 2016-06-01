% code for ECE HW1
% this script simulates the dynamics of a system of 3 metabolites (M)
% with known initial conditions and reaction dynamics controlled by
% 2 parameters (theta). system dynamics described in fcn 'system_dyn.m'
% Goals: (1) simulate and visualize dynamics
%        (2) simulate measurement values with noise
%        (3) use Bayesian inference to infer value of unknown parameter
%            from noisy data

%% (1) Simulate Dynamics

% initial concentrations in system (known)
M0 = [1; 1; 1];

% parameters: theta(2) = 2 (known). theta(1) is in [0,2] (unknown).
% here we assume the "true" value theta(1) = 0.9
theta = [0.9; 2];

% time range of simulation
trange = [0, 10];  % for solution at each solver timestep (smoother)
trange_m = 0:10;   % for solution only at specified sampling timepoints

% simulate with dynamics described in 'system_dyn.m'
% by solving initial value problem with ode45
[t, M] = ode45(@(t,M,theta)system_dyn(t,M,theta),trange,M0,[],theta);
[t_m, M_m] = ode45(@(t,M,theta)system_dyn(t,M,theta),trange_m,M0,[],theta);

% parametrically plot the concentration dynamics of all 3 metabolites
figure(1);
plot3(M(:,1),M(:,2),M(:,3)); hold all;
plot3(M_m(:,1),M_m(:,2),M_m(:,3));
xlabel('M1 concentration'); 
ylabel('M2 concentration'); 
zlabel('M3 concentration');
title('simulated system dynamics');
legend('smooth dynamics','discretely sampled dynamics','Location','NorthEast');
hold off;

%% (2) Simulate Measurement Values with Noise

% here we use the same values for M0 and theta and simulate again, making
% two changes: adding noise and measuring only at "discrete" time points.
% since the originally simulated system is deterministic, we can simply use
% the solutions for M3 above and add noise post-hoc.

% Gaussian noise parameters
n_mean = 0;
n_std = 1;

% create simulated noisy samples (M3hat)
M3hat = M_m(:,3) + randn(size(M_m(:,3))).*n_std + n_mean;
M3hat(1) = 0;        % assume no noise in known initial condition

% plot simulated true dynamics and noisy samples
figure(2); 
plot(t, M(:,3)); hold all;
plot(t_m, M3hat(:,1));
xlabel('time'); 
ylabel('M3 concentration'); 
title('dynamics of metabolite 3 (M3)');
legend('smooth simulated dynamics','noisy sampled dynamics','Location','NorthEast');
hold off;

%% (3) Infer "True" Theta(2) Value from Noisy Measurements

% we use Bayesian inference to calculate the posterior probability
% distribution over possible theta(1) values (in the range [0,2]).
% set of range of thetas to test (for theta(1); theta(2) is known):
n_theta = 21;        % number of test theta values
theta_test = [linspace(0,2,21); repmat(theta(2),1,n_theta)];

% we define a matrix to hold our prior (first row) and sequentially 
% calculated posteriors (subsequent rows) through a defined number of samples.
n_samp = 8;
post_mat = zeros(n_samp+1, n_theta);
post_mat(1,:) = 1/2; % uniform prior

% We now compare the dynamics of each system defined by the candidate
% test value for theta
for j=1:n_theta
    % For each index j, we take the j-th theta set to test and simulate its
    % dynamics, storing it in a temporary matrix Mtest
    [~,Mtest] = ode45(@(t,M,theta)system_dyn(t,M,theta),t_m,M0,[],theta_test(:,j));
    % Now, for each sequential time sample, we can update the posterior for 
    % theta by incorporating M3hat and calculating the likelihood wrt Mtest
    % NOTE that we do not calculate the normalizing evidence explicitly; we
    % instead treat each subsequent posterior as proportional to the
    % preceding one, and we will normalize all pdfs to 1 afterwards.
    for ts = 1:n_samp
        post_mat(ts+1,j)=post_mat(ts,j)*normpdf(M3hat(ts+1),Mtest(ts+1,3),n_std^2)/n_theta;
    end
end
% We now need to normalize the posteriors,
% also taking grid size into account
post_mat = bsxfun(@rdivide, post_mat, sum(post_mat)./n_theta);
post_mat(1,:) = 1/2; % uniform prior
%posterior_t2=posterior_t2/(sum(posterior_t2)/nsteps);
%posterior_t3=posterior_t3/(sum(posterior_t3)/nsteps);
%posterior_t4=posterior_t4/(sum(posterior_t4)/nsteps);

% plot the prior and subsequent posteriors after each time sample
figure(3);
plot(theta_test(1,:), post_mat); hold all;
title('sequentially updated posteriors for theta(2)');
xlabel('theta(2) value');
ylabel('posterior distribution');
legend('prior','posteriors');
hold off;

figure(4);
for ts=1:n_samp+1
    subplot(2,5,ts);
    plot(theta_test(1,:),post_mat(ts,:)); hold all;
    suptitle('sequentially updated posteriors for theta(2)');
    xlabel('theta(2) value');
    ylabel('distribution');
end
hold off;
