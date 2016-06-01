% group v. individual

%% GROUP PREDICTION
shelterSpace = 50;
groupSize = 25;
groupProtocol = 'mR';
N = 50;
p = 100;  % number of times to simulate
[~,~,~,~,missing] = load_evac_data(0);
toPlot = 1;

% train model on individual trials (output 'params')
trials = trial_ix('ind',shelterSpace,1,missing); bins = -0.05:0.1:1.05; 
[H,J,~,~,~,Q1,T1,P1,C1] = load_evac_data(0,trials,bins,ss); nts = size(Q1,2);
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
grp_dtime = zeros(numel(trials),N/groupSize,p);
grp_grpIDs = zeros(numel(trials),p,N);
grp_dtime_var = zeros(numel(trials),N/groupSize,p);
grp_dtime_dist_samp = zeros(numel(trials),N/groupSize,p);
grp_dtime_dist_empi = zeros(numel(trials),N/groupSize,p);
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
    if toPlot
    figure(1); hold on; 
    if groupSize==5
     if strcmp(groupProtocol,'fTG')
      if tr==trials(end)
        subplot(3,4,[7,8,11,12]);
      end
      if find(trials==tr)<7 
        subplot(3,4,find(trials==tr));
      end
      if find(trials==tr)==7 || find(trials==tr)==8
        subplot(3,4,find(trials==tr)+2);
      end
     end
     if strcmp(groupProtocol,'lTG')
      if tr==trials(end)
        subplot(2,4,[3,4,7,8]);
      end
      if find(trials==tr)<3 
        subplot(2,4,find(trials==tr));
      end
      if find(trials==tr)==3 || find(trials==tr)==4
        subplot(2,4,find(trials==tr)+2);
      end
     end
    end
    if groupSize==25
     if strcmp(groupProtocol,'mR')
      if tr==trials(end)
        subplot(3,3,[7,8,11,12]);
      end
      if find(trials==tr)<7 
        subplot(3,3,find(trials==tr));
      end
      if find(trials==tr)==7 || find(trials==tr)==8
        subplot(3,3,find(trials==tr)+2);
      end
     end
     if strcmp(groupProtocol,'lTG')
      if tr==trials(end)
        subplot(2,4,[3,4,7,8]);
      end
      if find(trials==tr)<3 
        subplot(2,4,find(trials==tr));
      end
      if find(trials==tr)==3 || find(trials==tr)==4
        subplot(2,4,find(trials==tr)+2);
      end
     end
    end
            plot(C1(tr,:),'k--'); hold on; % observed in experiment
            plot(Q1(tr,:).*N,':');         % Phit trajectory
            for i=1:10%p
                plot(tgrp(:,i),Cgrp(:,i),'-o'); hold all;
                %plot(t_mat(r_mat(:,i)~=0,i),ix(r_mat(r_mat(:,i)~=0,i))-1,':o'); hold all;
            end
            plot(tbins,mean(Cbin,2),'-k');
            title(['trial ' num2str(tr) ', ' num2str(p) ' samples']); axis([0 60 0 50])
            if tr==trials(end)
            legend('empirical data','Phit trajectory (scaled by N)',...
                  [groupProtocol ' sampled trajectories'],...
                   'location','northwest');
            end
            xlabel('time'); ylabel('cumulative no. evacuated');
          hold off;
    end
    
    % for each relevant trial, compute nevac_b and nevac(t)
    nevac_t(trials==tr,1:trEnd(tr)) = mean(Cbin(2:trEnd(tr)+1,:),2)' - C1(tr,1:trEnd(tr));
    nevac_b(trials==tr) = mean(Cbin(trEnd(tr)+1,:)) - C1(tr,end);            
    
    % make a plot that breaks down group dynamics in a single instance
    if toPlot
    px = 1;
    t_sort = sortentry([t_evac(:,px);NaN*zeros(N-size(t_evac,1),1)],'col',0,grpIDs(px,:));
    t_sort = reshape(t_sort,groupSize,numel(unique(grpIDs)));
    switch groupProtocol
        case 'fTG',
          t_sort = sortentry(t_sort,'row',0,t_sort(1,:));
        case 'lTG',
          t_sort = sortentry(t_sort,'row',0,t_sort(end,:));
        case 'mR',
          t_sort = sortentry(t_sort,'row',0,t_sort(floor(groupSize/2)+1,:));  
        case 'ind',
          t_sort = sortentry(t_sort,'row',0,t_sort(1,:));
    end
    figure(2); hold on;
        if strcmp(groupProtocol,'fTG')
          if tr==trials(end)
            subplot(3,4,[7,8,11,12]);
          end
          if find(trials==tr)<7 
            subplot(3,4,find(trials==tr));
          end
          if find(trials==tr)==7 || find(trials==tr)==8
            subplot(3,4,find(trials==tr)+2);
          end
        end
        if strcmp(groupProtocol,'lTG')
          if tr==trials(end)
            subplot(2,4,[3,4,7,8]);
          end
          if find(trials==tr)<3 
            subplot(2,4,find(trials==tr));
          end
          if find(trials==tr)==3 || find(trials==tr)==4
            subplot(2,4,find(trials==tr)+2);
          end
        end
            plot(C1(tr,:),'k--'); hold all;  % observed in experiment
            plot(Q1(tr,:).*N,':');         % Phit trajectory
            plot(tbins,mean(Cbin,2),'-k');
            for grp = 1:numel(unique(grpIDs))
                plot(t_sort(:,grp),(1:groupSize)+(grp-1)*groupSize,'-');
            end
            for grp = 1:numel(unique(grpIDs))
                plot(1:size(Q1,2),ones(size(Q1,2),1).*(grp-1)*groupSize,':k');
            end
            title(['trial ' num2str(tr) ', ' num2str(p) ' samples']); axis([0 60 0 50])
            if trials(end)==tr
            legend('empirical data','Phit trajectory (scaled by N)',...
                   'mean sampled evacuations',...
                   'individual decisions',...
                   'location','northwest');
            end
            xlabel('time'); ylabel('cumulative no. evacuated');
    end
    
    % for each group, in each instance,
    % compute std in evac time, and ind diffs from actual grp evac times            
	dtime_var = zeros(numel(unique(grpIDs)),p);
    dtime_dist_samp = zeros(numel(unique(grpIDs)),p);
    dtime_dist_empi = zeros(numel(unique(grpIDs)),p);
    empi_times = [];
        for i=find(diff(C1(tr,:)))
            empi_times = [empi_times;repmat(i,diff(C1(tr,i:i+1))/groupSize,1)];
        end
        try assert(length(empi_times)==numel(unique(grpIDs)));
        catch
        	empi_times = [empi_times;NaN*ones(numel(unique(grpIDs))-length(empi_times),1)];
        end
    for px=1:p
      % sort all individual evac times into groups  
      t_sort = sortentry([t_evac(:,px);NaN*zeros(N-size(t_evac,1),1)],'col',0,grpIDs(px,:));
      t_sort = reshape(t_sort,groupSize,numel(unique(grpIDs)));
      % now sort all group lists of evac times by time of whole grp evac
      switch groupProtocol
        case 'fTG',
          [t_sort,grp_ix] = sortentry(t_sort,'row',0,t_sort(1,:));
        case 'lTG',
          [t_sort,grp_ix] = sortentry(t_sort,'row',0,t_sort(end,:));
        case 'mR',
          [t_sort,grp_ix] = sortentry(t_sort,'row',0,t_sort(floor(groupSize/2)+1,:));  
        case 'ind',
          [t_sort,grp_ix] = sortentry(t_sort,'row',0,t_sort(1,:));
      end
      [~,new_ix] = sort(grp_ix);
      grpIDs_new = zeros(size(grpIDs));
      % calculate stats on ind v. group evac times
      for grp = 1:numel(unique(grpIDs))
        gtimes = t_sort(:,grp);
        % quantify variance in decision time
        dtime_var(grp,px) = nanstd(gtimes);
        % quantify mean distance from sampled group evac time
        dtime_dist_samp(grp,px) = sqrt(nanmean((tgrp(grp,px)-gtimes).^2));
        dtime_dist_empi(grp,px) = sqrt(nanmean((empi_times(grp)-gtimes).^2));
        % reorder grpIDs to reflect grp evac time order
        grpIDs_new(px,grpIDs(px,:)==grp) = new_ix(grp);
      end
    end  % end loop over instances (px)
    
% plot the variances and distances from simulated and actual evac times    
if toPlot
figure(3); hold on;
    if strcmp(groupProtocol,'fTG')
          if tr==trials(end)
            subplot(3,4,[7,8,11,12]);
          end
          if find(trials==tr)<7 
            subplot(3,4,find(trials==tr));
          end
          if find(trials==tr)==7 || find(trials==tr)==8
            subplot(3,4,find(trials==tr)+2);
          end
    end
    if strcmp(groupProtocol,'lTG')
          if tr==trials(end)
            subplot(2,4,[3,4,7,8]);
          end
          if find(trials==tr)<3 
            subplot(2,4,find(trials==tr));
          end
          if find(trials==tr)==3 || find(trials==tr)==4
            subplot(2,4,find(trials==tr)+2);
          end
    end
        hist(dtime_var');
        %xlabel('group evacuation time standard deviation');
        xlabel('std (time steps)');
        ylabel('number of groups');
        %suptitle(['Spread in Group Evacuation Times,' groupProtocol ', ss=' num2str(shelterSpace) ', gs=' num2str(groupSize)]);
        leg = {'grp1','grp2','grp3','grp4','grp5','grp6','grp7','grp8','grp9','grp10'};
        if trials(end)==tr
        legend(leg);
        end
        hold off;
% figure; hist(dtime_dist_samp');
%         xlabel('difference between ind. and group evac times (time steps)');
%         ylabel('number of groups');
%         title(['Sampled Individual v. Simulated Group Evac Times,' groupProtocol ', ss=' num2str(shelterSpace) ', gs=' num2str(groupSize)]);
%         legend(leg);
figure(4); hold on;
    if strcmp(groupProtocol,'fTG')
          if tr==trials(end)
            subplot(3,4,[7,8,11,12]);
          end
          if find(trials==tr)<7 
            subplot(3,4,find(trials==tr));
          end
          if find(trials==tr)==7 || find(trials==tr)==8
            subplot(3,4,find(trials==tr)+2);
          end
    end
    if strcmp(groupProtocol,'lTG')
          if tr==trials(end)
            subplot(2,4,[3,4,7,8]);
          end
          if find(trials==tr)<3 
            subplot(2,4,find(trials==tr));
          end
          if find(trials==tr)==3 || find(trials==tr)==4
            subplot(2,4,find(trials==tr)+2);
          end
    end
        hist(dtime_dist_empi');
        xlabel('RMS difference (time steps)');
        %xlabel('difference between ind. and group evac times (time steps)');
        ylabel('number of groups');
        %suptitle(['Sampled Individual v. Empirical Group Evac Times,' groupProtocol ', ss=' num2str(shelterSpace) ', gs=' num2str(groupSize)]);
        if trials(end)==tr
        legend(leg);
        end
        hold off;
end

grp_grpIDs(trials==tr,:,:) = shiftdim(grpIDs_new,-1);
grp_dtime(trials==tr,:,:) = shiftdim(tgrp,-1);
grp_dtime_var(trials==tr,:,:) = shiftdim(dtime_var,-1);
grp_dtime_dist_samp(trials==tr,:,:) = shiftdim(dtime_dist_samp,-1);
grp_dtime_dist_empi(trials==tr,:,:) = shiftdim(dtime_dist_empi,-1);
end  % end loop over trials
figure(3); hold on;
suptitle(['Spread in Group Evacuation Times,' groupProtocol ', ss=' num2str(shelterSpace) ', gs=' num2str(groupSize)]);
hold off;
figure(4); hold on;
suptitle(['Sampled Individual v. Empirical Group Evac Times,' groupProtocol ', ss=' num2str(shelterSpace) ', gs=' num2str(groupSize)]);
hold off;

%figure(1); hold on; tightfig; hold off;
%figure(2); hold on; tightfig; hold off;
%figure(3); hold on; tightfig; hold off;
%figure(4); hold on; tightfig; hold off;

%% win/lose strategies: plot stats
% info: did each group evacuate?
didEvacs = (grp_dtime<=repmat(trEnd(trials)',1,N/groupSize,p));
% info: did each trial hit?
didHits = Q1(trials,end);

% stats by hit/miss
varsHit = grp_dtime_var(~~didHits,:,:);
varsMiss = grp_dtime_var(~didHits,:,:);
distHit = grp_dtime_dist_empi(~~didHits,:,:);
distMiss = grp_dtime_dist_empi(~didHits,:,:);
nbins = 15;
figure; bar_unpaired(2,varsHit,varsMiss,'nbins',nbins);
        legend(['hit trials (' num2str(sum(didHits)) ')'],['miss trials (' num2str(sum(~didHits)) ')']);
        xlabel('variance in group evac time');
        ylabel('number of groups');
        title(['grp evac time variances, ' groupProtocol ', ss=' num2str(shelterSpace) ', groupSize=' num2str(groupSize)]);
figure; bar_unpaired(2,distHit,distMiss,'nbins',nbins);
        legend(['hit trials (' num2str(sum(didHits)) ')'],['miss trials (' num2str(sum(~didHits)) ')']);
        xlabel('RMS difference between ind. simulated and group empirical evac times');
        ylabel('number of groups');
        title(['ind. sim v. grp emp evac time differences, ' groupProtocol ', ss=' num2str(shelterSpace) ', groupSize=' num2str(groupSize)]);

% plot stats as a function of group success
varsHE = varsHit(~~didEvacs(~~didHits,:,:));
varsME = varsMiss(~~didEvacs(~didHits,:,:));
varsHN = varsHit(~didEvacs(~~didHits,:,:));
varsMN = varsMiss(~didEvacs(~didHits,:,:));
figure; bar_unpaired(4,varsHE,varsMN,varsHN,varsME,'nbins',nbins);
        legend(['hit (' num2str(sum(didHits)) '), evac'],['miss (' num2str(sum(~didHits)) '), no evac'],'hit, no evac','miss, evac');
        xlabel('variance in group evac time');
        ylabel('number of groups');
        title(['grp evac time variances, ' groupProtocol ', ss=' num2str(shelterSpace) ', groupSize=' num2str(groupSize)]);
distHE = distHit(~~didEvacs(~~didHits,:,:));
distME = distMiss(~~didEvacs(~didHits,:,:));
distHN = distHit(~didEvacs(~~didHits,:,:));
distMN = distMiss(~didEvacs(~didHits,:,:));
figure; bar_unpaired(4,distHE,distMN,distHN,distME,'nbins',nbins);
        legend(['hit (' num2str(sum(didHits)) '), evac'],['miss (' num2str(sum(~didHits)) '), no evac'],'hit, no evac','miss, evac');
        xlabel('RMS difference between ind. simulated and group empirical evac times');
        ylabel('number of groups');
        title(['ind. sim v. grp emp evac time differences, ' groupProtocol ', ss=' num2str(shelterSpace) ', groupSize=' num2str(groupSize)]);

        
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

% for each group on each trial...
grpScore = zeros(ntrials);

% is there a correlation between GROUP rank and....
% (a) mean individual rank within group

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        