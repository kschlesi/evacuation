function [mean_evac,T] = mastereq(decmodel,tf)
tic;
%% initialize generator matrix

% total number of states =  # individuals + 1
N = 51;

% initial probability distribution 
% all individuals are at home
P_init = zeros(N,1);
P_init(1) = 1;

% each run is 60 seconds
% tf = 60;
T_range = [0 tf];
% define an evenly-spaced grid of time points, where we have data at each
% time point
qt=linspace(0,tf,61);

% q is the probability of evacuation as a function of the disaster
% probability
q = decmodel;

% A is the generator matrix that multiplies the probability vector in the
% master equation
% we determine A at every time step
A = zeros(N,N,length(qt)); 

for k=1:length(qt)
    for n=1:N
        for i=1:N
            if i < n
                A(n,i,k)=nchoosek(N-i,n-i)*q(k)^(n-i)*(1-q(k))^(N-n);
            end
        end
    end
    for i=1:N
        A(i,i,k)=-sum(A(:,i,k)); % the diagonal entries of A are the negative sums of the respective columns
    end
end
Aa=reshape(A,[],length(q)).'; 

%% integrate master equation
[T, P] = ode45(@(t,P) odefunc(P,qt,Aa,t,N), T_range, P_init, []);

% calculate mean number evacuated
mean_evac = P*(0:(N-1))';
toc;
%% 
function dP=odefunc(P,qt,Aa,t,N)
A_interp=interp1(qt,Aa,t); % interpolate values of generator matrix at ODE-solver timesteps
finalA=reshape(A_interp,N,N); 
dP=finalA*P;
end

end