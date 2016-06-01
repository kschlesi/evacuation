% GILLESPIE ALGORITHM IMPLEMENTATION

T = 10;       % duration of reaction to simulate
t0 = 0;       % time to start
nS = 100;     % number of simulations to run
tstep = 0.01;  % time interval at which to sample results for averaging

n0 = [10;     % initial count of C1
       5];    % initial count of C2

theta = [10;  % reaction rates
          1;
         50;
          1];

nu = [1,  1, -1, -1;  % state-change matrix
      0, -1,  1,  0];
  
a = @(n,theta)[ theta(1);        % reaction propensity functions
                theta(2)*n(1)*(n(1)-1)*n(2)/2;
                theta(3)*n(1);
                theta(4)*n(1) ];

% run simulations
[t_mat,n_mat,r_mat,ns_mat] = gillespie(theta,a,nu,n0,t0,T,nS);

% sample from consistent timepoints for each simulation
t_samp = t0:tstep:t0+T;  % time points at which to sample
tx_mat = zeros(length(t_samp),nS);  % sampled indices in each simulation
n_samp = zeros(length(n0),length(t_samp),nS);
r_samp = zeros(length(t_samp),nS);
for sim=1:nS
  for tt = t_samp
    tx_mat(t_samp==tt,sim) = find(t_mat(1:ns_mat(sim),sim)<=tt,1,'last');
  end
  n_samp(:,:,sim) = n_mat(:,tx_mat(:,sim),sim);
  r_samp(:,sim) = r_mat(tx_mat(:,sim),sim);
end


% figure: plot time course for simulations 1-3
nR = length(theta);
for sim=1:min(nS,3)
    rge = 1:ns_mat(sim);
    figure;
    [ax,~,p2] = plotyy(t_mat(rge,sim),n_mat(:,rge,sim),...
           t_mat(rge,sim),r_mat(rge,sim).*max(max(max(n_mat)))./nR);
    hold on; p2.Marker = 'o'; p2.LineStyle = 'none';
    ax(2).YTick = linspace(min(ax(2).YTick),max(max(max(n_mat))),nR+1);
    ax(2).YTickLabel = str2double(ax(2).YTickLabel)./max(max(max(n_mat))).*nR;
    ylabel(ax(1),'population count');
    ylabel(ax(2),'reaction ID');
    xlabel('time');
    legend('C1','C2','rxn');

%     figure;
%     plot(t_mat(rge,sim)); hold all;
%     plot(r_mat(rge,sim).*max(t_mat(rge,sim))./4,'o');
%     xlabel('reaction count');
%     ylabel('time');
end

% figure: plot sampled time courses for sims 1-3
for sim=1:min(nS,3)
figure;
    [ax,~,p2] = plotyy(t_samp,n_samp(:,:,sim),...
                        t_samp,r_samp(:,sim).*max(max(max(n_samp)))./nR);
    hold on; p2.Marker = 'o'; p2.LineStyle = 'none';
    ax(2).YTick = linspace(min(ax(2).YTick),max(max(max(n_samp))),nR+1);
    ax(2).YTickLabel = str2double(ax(2).YTickLabel)./max(max(max(n_samp))).*nR;
    ylabel(ax(1),'population count');
    ylabel(ax(2),'reaction ID');
    xlabel('time');
    legend('C1','C2','rxn');
end

% figure: plot average of all sampled time courses
figure;
    [ax,~,p2] = plotyy(t_samp,mean(n_samp,3),...
                        t_samp,mode(r_samp,2).*max(max(max(n_samp)))./nR);
    hold on; p2.Marker = 'o'; p2.LineStyle = 'none';
    ax(2).YTick = linspace(min(ax(2).YTick),max(max(max(n_samp))),nR+1);
    ax(2).YTickLabel = str2double(ax(2).YTickLabel)./max(max(max(n_samp))).*nR;
    ylabel(ax(1),'population count');
    ylabel(ax(2),'reaction ID');
    xlabel('time');
    legend('C1','C2','rxn');        
    
% figure: switch from regime 1 to 2
n1ax = linspace(1,10,1000);
n2ax = 2*(theta(1)+(theta(3)+theta(4)).*n1ax)./(n1ax.*(n1ax-1));
n1ax0 = 1:1:10;
n2ax0 = 2*(theta(1)+(theta(3)+theta(4)).*n1ax0)./(n1ax0.*(n1ax0-1));
figure; plot(n1ax,n2ax); hold on; plot(n1ax0,n2ax0,'o');
axis([1 10 0 125]);
xlabel('n_1'); ylabel('n_2');
legend('a_2 = a_1+a_3+a_4');


