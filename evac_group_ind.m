% group v. individual

%% GROUP PREDICTION
shelterSpace = 50;
groupSize = 5;
groupProtocol = 'fTG';
N = 50;
p = 100;  % number of times to simulate
[~,~,~,~,missing] = load_evac_data(0);

% train model on individual trials (output 'params')
trials = trial_ix('ind',shelterSpace,1,missing); bins = -0.05:0.1:1.05; 
[H,J,~,~,~,Q1,T1,P1,C1] = load_evac_data(0,trials,bins); nts = size(Q1,2);
% qform = @(Phit_,pv_) pv_(1).*(Phit_.^pv_(3))./(pv_(2).^pv_(3)+Phit_.^pv_(3));
% Phit = 0:0.1:1;
% qfit = @(pv_) qform(Phit(:,2:end),pv_);
% MLfit = @(pvec) -1*sum(((sum(H(:,2:end))-sum(J(:,2:end))).*log(1-qfit(pvec)) + sum(J(:,2:end)).*log(qfit(pvec))));
% startp = [1;0.5;10];
qform = @(Phit_,pv_) pv_(1).*Phit_.^pv_(2);
Phit = 0:0.1:1;
qfit = @(pv_) qform(Phit(:,2:end-1),pv_);
MLfit = @(pvec_) -1*sum(((sum(H(:,2:end-1))-sum(J(:,2:end-1))).*log(1-qfit(pvec_)) + sum(J(:,2:end-1)).*log(qfit(pvec_))));
startp = [1.1;1];
options = optimoptions(@fminunc,'MaxFunEvals',10000);
[params1,MLval] = fminunc(MLfit,startp,options);
A = zeros(length(startp)); b = zeros(1,length(startp)); q_con = @(pvec_)(qform(Phit,pvec_)-1); 
[theta_con,MLval_con] = fmincon(MLfit,startp,A,b,A,b,0,100,@(pvec_)q_add_eq(pvec_,q_con));
if MLval<MLval_con && all(qform(Phit,params1)<=1)
    params = params1;
else
    params = theta_con;
end

% choose trials to predict
trials = trial_ix(groupProtocol,shelterSpace,groupSize,missing);
plotH = trials(find(Q1(trials,end)==1,1,'first'));
plotM = trials(find(Q1(trials,end)==0,1,'first'));
trEnd = zeros(1,size(Q1,1)); % we will store the time at which the trials ended
for j = 1:size(Q1,1) % iterate through each trial
    trEnd(j) = find(Q1(j,:)==Q1(j,end),1,'first');
end
nevac_b = zeros(numel(trials),1);
nevac_t = zeros(numel(trials),nts);
for tr = trials(:)'
    % randomly generate N individual evac times, p simulations
    % for comments see similar section above
    q = @(Phit_) qform(Phit_,params);
    tbins = 0:1:60;
    a = @(t_,P_,theta_) a_vec(t_,P_,theta_{1},theta_{2}); % this returns a prop. vector of only lower-triangular vals
    aix = tril(ones(N+1),-1);  % indices of lower-triangular entries
    theta = {makeA(N,q(Phit)); Q1(tr,:)}; % these are the non-explicitly-time-dependent input parameters
    P0 = [1;zeros(N,1)];   % this is the initial condition
    t0 = 0;                % starting time of simulation
    endx = find(Q1(tr,:)==Q1(tr,end),1,'first');
    T = endx + 1 - (endx==length(Q1(tr,:)));  % duration of simulation
    [ix,jx] = ind2sub(size(aix),find(~~aix));
    sx = sub2ind(size(aix),find(~~aix));
    nu = sparse(jx,1:length(sx),-1,N+1,length(sx)) + sparse(ix,1:length(sx),1);
    
    % this generates the samples with Gillespie algorithm...
    [t_mat,n_mat,r_mat,ns_mat,br_mat] = gillespie_tdep(theta,a,nu,P0,t0,T,tbins,p);
    
    % this creates a list of evac times & group IDs
    % randomly generate groups of size groupSize
    grpID1 = repmat(1:floor(N/groupSize),1,groupSize);
    grpID1 = [grpID1 zeros(1,N-length(grpID1))];
    grpIDs = zeros(p,N);
    n_evac = zeros(size(r_mat)); 
    n_evac(~~r_mat) = ix(r_mat(~~r_mat))-jx(r_mat(~~r_mat));
    t_evac = zeros(max(sum(n_evac)),p).*NaN;
    for px=1:p
      tx = [];
      for ne = removeval(unique(n_evac),0)'
         tx = [tx; repmat(t_mat(n_evac(:,px)==ne,px),ne,1)];
      end
      t_evac(1:sum(n_evac(:,px)),px) = sort(tx)';
      grpIDs(px,:) = grpID1(randperm(N));
    end
    
    % this creates cumulative evac plots given group IDs and evac times
    [Cgrp,tgrp,Cbin] = cum_evac(t_evac,grpIDs,groupProtocol,tbins);
    % plot cumulative evacuated
    figure; plot(C1(tr,:),'k--'); hold on; % observed in experiment
            plot(Q1(tr,:).*N,':');         % Phit trajectory
            for i=1:10%p
                plot(tgrp(:,i),Cgrp(:,i),'-o'); hold all;
                %plot(t_mat(r_mat(:,i)~=0,i),ix(r_mat(r_mat(:,i)~=0,i))-1,':o'); hold all;
            end
            plot(tbins,mean(Cbin,2),'-k');
            title(['trial ' num2str(tr) ', ' num2str(p) ' samples']); axis([0 60 0 50])
            legend('empirical data','Phit trajectory (scaled by N)',...
                   groupProtocol,...%'individual',
                   'location','northwest');
            xlabel('time'); ylabel('cumulative no. evacuated');
    
    % for each relevant trial, compute nevac_b and nevac(t)
    nevac_t(trials==tr,1:trEnd(tr)) = mean(Cbin(2:trEnd(tr)+1,:),2)' - C1(tr,1:trEnd(tr));
    nevac_b(trials==tr) = mean(Cbin(trEnd(tr)+1,:)) - C1(tr,end);            
    
    % make a plot that breaks down group dynamics in a single instance
    px = 1;
    t_sort = sortentry([t_evac(:,px);NaN*zeros(N-size(t_evac,1),1)],'col',0,grpIDs(px,:));
    figure; plot(C1(tr,:),'k--'); hold on; % observed in experiment
            plot(Q1(tr,:).*N,'.--');         % Phit trajectory
            plot(tbins,mean(Cbin,2),'-k');
            for grp = 1:numel(unique(grpIDs))
                plot(t_sort((grp-1)*groupSize+1:grp*groupSize),(1:groupSize)+(grp-1)*groupSize,'-');
            end
            for grp = 1:numel(unique(grpIDs))
                plot(1:size(Q1,2),ones(size(Q1,2),1).*(grp-1)*groupSize,':k');
            end
            title(['trial ' num2str(tr) ', ' num2str(p) ' samples']); axis([0 60 0 50])
            legend('empirical data','Phit trajectory (scaled by N)',...
                   'mean sampled evacuations',...
                   'individual decisions',...
                   'location','northwest');
            xlabel('time'); ylabel('cumulative no. evacuated');
    
end  % end loop over trials

%% plots of overall trialwise error stats and performance

% plot nevac_b for each game
horm = Q1(:,end);
hits = NaN*zeros(numel(trials),1); misses = hits;
hits(~~horm(trials)) = nevac_b(~~horm(trials)); 
misses(~horm(trials)) = nevac_b(~horm(trials));
figure; plot(nevac_b); hold on; 
        plot(zeros(numel(trials),1),'--k');
        plot(hits,'or'); plot(misses,'ok');
        xlabel('trial number'); set(gca,'XTickLabel',trials);
        ylabel('error in total number evacuated');
        title(['number evacuated prediction error, ' groupProtocol ', ss=' ...
               num2str(shelterSpace) ', grp=' num2str(groupSize)]);
           
% save from several groupProtocols with...
%nevac_t_fTG_50_5 = nevac_t; 
%nevac_b_fTG_50_5 = nevac_b;
eval(['nevac_t_' groupProtocol '_' num2str(shelterSpace) '_' num2str(groupSize) ' = nevac_t;']);
eval(['nevac_b_' groupProtocol '_' num2str(shelterSpace) '_' num2str(groupSize) ' = nevac_tb;']);

% plot & compare
figure;
for gP = 3
    switch gP
        case 1, toPlot = nevac_t_fTG_50_5;
                mean_abs_b(1,1) = mean(abs(nevac_b_fTG_50_5));
                mean_abs_b(2,1) = mean(mean(abs(nevac_t_fTG_50_5)));
                plot(1:60,toPlot',':b'); hold on;
                plot(1:60,mean(toPlot),'-b');
        case 2, toPlot = nevac_t_lTG_50_5;
                mean_abs_b(1,2) = mean(abs(nevac_b_lTG_50_5));
                mean_abs_b(2,2) = mean(mean(abs(nevac_t_lTG_50_5)));
                plot(1:60,toPlot',':m'); hold on;
                plot(1:60,mean(toPlot),'-m');
        case 3, toPlot = nevac_t_mR_50_5;
                mean_abs_b(1,3) = mean(abs(nevac_b_mR_50_5));
                mean_abs_b(2,3) = mean(mean(abs(nevac_t_mR_50_5)));
                plot(1:60,toPlot'); hold on;
                plot(1:60,toPlot'); hold on;
                plot(1:60,toPlot,'-');
    end
end
title(['number evacuated prediction error, ss=' ...
               num2str(shelterSpace) ', grp=' num2str(groupSize)]);
xlabel('time'); ylabel('number evacuated error');
%legend('fTG','lTG','mR');

figure; bar(mean_abs_b');
        xlabel('group protocol'); set(gca,'XTickLabel',{'fTG';'lTG';'mR'});
        ylabel('mean absolute error');
        title(['total evacuations, mean absolute prediction error, ss=' ...
               num2str(shelterSpace) ', grp=' num2str(groupSize)]);
        legend('trial end error','mean over timesteps');
        
%% looking at the DATA itself...

% is there a correlation between GROUP rank and....
% (a) 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        