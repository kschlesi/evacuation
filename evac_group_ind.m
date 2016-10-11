% group v. individual

%% GROUP PREDICTION
shelterSpace = 25;
groupSize = 5;
groupProtocol = 'fTG';
N = 50;
p = 100;  % number of times to simulate
[~,~,~,~,missing] = load_evac_data(0);
toPlot = false;
toPlot_examples = true;

switch groupProtocol
                case 'fTG'
                  gpstr = 'FTG';
                case 'mR'
                  gpstr = 'MV';
                case 'lTG'
                  gpstr = 'LTG';
end

% train model on individual trials (output 'params')
trials = trial_ix('ind',shelterSpace,1,missing); bins = 0:0.1:1.1; %bins = -0.05:0.1:1.05; 
[H,J,~,~,~,Q1,T1,P1,C1] = load_evac_data(0,trials,bins); nts = size(Q1,2);

qform = @(Phit_,pv_) pv_(1).*Phit_.^pv_(2);
startp = [1.1;1];
Phit = 0:0.1:1;
params = ML_fit_beta(qform,Phit(2:end-1),H(:,2:end-1),J(:,2:end-1),startp);

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
    %%% fix As for ss<N
    ss = shelterSpace;
    if ss<N
      for i=1:length(Phit)
        % zero blocked transitions and transfer probabilities to allowed ones
        theta{1}(ss+1,:,i) = theta{1}(ss+1,:,i) + sum(theta{1}(ss+2:end,:,i));
        theta{1}(ss+2:end,:,i) = 0; 
        theta{1}(:,ss+2:end,i) = 0;
        % zero and re-compute diagonal
        theta{1}(:,:,i) = theta{1}(:,:,i).*(~eye(size(theta{1}(:,:,i)))); assert(all(~diag(theta{1}(:,:,i))));
        theta{1}(:,:,i) = theta{1}(:,:,i) - diag(sum(theta{1}(:,:,i))); %assert(all(~sum(A)));
      end
    end
        
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
    [Cgrp,tgrp,Cbin] = cum_evac(t_evac,grpIDs,groupProtocol,tbins,'ss',shelterSpace);
    % plot cumulative evacuated
    if toPlot_examples
    figure(1); hold all;
          if groupSize==5 && shelterSpace==50
            makeLeg = false;
            switch groupProtocol
                case 'fTG',
                  subplot(3,6,find(trials==tr));
                case 'mR'
                  subplot(3,6,[11,12,17,18]);
                  makeLeg = true;
                case 'lTG'
                  if find(trials==tr)<2  
                    subplot(3,6,find(trials==tr)+9);
                  else
                    subplot(3,6,find(trials==tr)+11);  
                  end
            end
          end
          if groupSize==25 && shelterSpace==50
             makeLeg = false;
             switch groupProtocol
                case 'fTG',
                  subplot(4,5,find(trials==tr));
                case 'mR'
                  if find(trials==tr)==1  
                  subplot(4,5,5);
                  else
                      if find(trials==tr)<7
                          subplot(4,5,find(trials==tr)+4);
                      else
                          subplot(4,5,[14,15,19,20]);
                          makeLeg = true;
                      end
                  end
                case 'lTG'
                  if find(trials==tr)<4
                    subplot(4,5,find(trials==tr)+10);
                  else
                    subplot(4,5,find(trials==tr)+12);
                  end
            end  
          end
          if groupSize==5 && shelterSpace==25
            makeLeg = false;
            switch groupProtocol
                case 'fTG',
                  subplot(3,7,find(trials==tr));
                case 'mR'
                  if find(trials==tr)<6
                    subplot(3,7,find(trials==tr)+7);
                  else
                    makeLeg = true;
                    subplot(3,7,[13,14,20,21]);
                  end
                case 'lTG'
                  subplot(3,7,find(trials==tr)+14);
            end
          end
          if groupSize==25 && shelterSpace==25
             makeLeg = false;
             switch groupProtocol
                case 'fTG',
                  subplot(4,5,find(trials==tr));
                case 'mR'
                  if find(trials==tr)<3  
                    subplot(4,5,find(trials==tr)+3);
                  else
                    subplot(4,5,find(trials==tr)+3);  
                  end
                case 'lTG'
                  if find(trials==tr)<4  
                    subplot(4,5,find(trials==tr)+10);
                  else
                    if find(trials==tr)<7  
                      subplot(4,5,find(trials==tr)+12);
                    else
                      makeLeg = true;
                      subplot(4,5,[14,15,19,20]);  
                    end
                  end
            end  
          end
    	    plot(C1(tr,:),'--','color',[0 0.5 0],'LineWidth',2);  hold all; % Phit trajectory
            plot((floor(Q1(tr,:)*10)/10).*N,'--','LineWidth',1);  % observed in experiment
	        plot(tbins,mean(Cbin,2),'-k','LineWidth',2); 
            for i=1:10%p
                plot(tgrp(:,i),Cgrp(:,i),':','LineWidth',1); 
                %plot(t_mat(r_mat(:,i)~=0,i),ix(r_mat(r_mat(:,i)~=0,i))-1,':o'); hold all;
            end
            %ax1.YColor = 'b';
            title(['trial ' num2str(tr) ', ' gpstr]); axis([0 60 0 50]);
            if makeLeg
                legend('empirical evacuations','Phit trajectory (scaled by N)',...
                   ['simulated evacs (mean over ' num2str(p) ' samples)'] ,...
                   ['randomly selected ' gpstr ' samples'],...%'individual',
                   'location','northwest');
            end
            xlabel('time step'); ylabel('cumulative no. evacuated');
            title(trial_conv(tr,ss,groupSize,groupProtocol));
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
    figure; plot(C1(tr,:),'k--'); hold all;  % observed in experiment
            plot(Q1(tr,:).*N,'.--');         % Phit trajectory
            plot(tbins,mean(Cbin,2),'-k');
            for grp = 1:numel(unique(grpIDs))
                plot(t_sort(:,grp),(1:groupSize)+(grp-1)*groupSize,'-');
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
figure; hist(dtime_var');
        xlabel('group evacuation time standard deviation');
        ylabel('number of groups');
        title(['Spread in Group Evacuation Times,' groupProtocol ', ss=' num2str(shelterSpace) ', gs=' num2str(groupSize)]);
        leg = {'grp1','grp2','grp3','grp4','grp5','grp6','grp7','grp8','grp9','grp10'};
        legend(leg);
% figure; hist(dtime_dist_samp');
%         xlabel('difference between ind. and group evac times (time steps)');
%         ylabel('number of groups');
%         title(['Sampled Individual v. Simulated Group Evac Times,' groupProtocol ', ss=' num2str(shelterSpace) ', gs=' num2str(groupSize)]);
%         legend(leg);
figure; hist(dtime_dist_empi');
        xlabel('difference between ind. and group evac times (time steps)');
        ylabel('number of groups');
        title(['Sampled Individual v. Empirical Group Evac Times,' groupProtocol ', ss=' num2str(shelterSpace) ', gs=' num2str(groupSize)]);
        legend(leg);
end

grp_grpIDs(trials==tr,:,:) = shiftdim(grpIDs_new,-1);
grp_dtime(trials==tr,:,:) = shiftdim(tgrp,-1);
grp_dtime_var(trials==tr,:,:) = shiftdim(dtime_var,-1);
grp_dtime_dist_samp(trials==tr,:,:) = shiftdim(dtime_dist_samp,-1);
grp_dtime_dist_empi(trials==tr,:,:) = shiftdim(dtime_dist_empi,-1);
end  % end loop over trials

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
grpRank = zeros(ntrials);

% is there a correlation between GROUP rank and....
% (a) mean individual rank within group

        
% group performance?


        
        
        
        if groupSize==5 && shelterSpace==50
            switch groupProtocol
                case 'fTG',
                  subplot(3,5,find(trials==tr));
                case 'mR'
                  subplot(3,5,10);
                case 'lTG'
                  subplot(3,5,find(trials==tr)+10);
            end
          end
          if groupSize==25 && shelterSpace==50
             switch groupProtocol
                case 'fTG',
                  subplot(4,4,find(trials==tr));
                case 'mR'
                  subplot(4,4,find(trials==tr)+4);
                case 'lTG'
                  subplot(4,4,find(trials==tr)+11);
            end  
          end
          if groupSize==5 && shelterSpace==25
            switch groupProtocol
                case 'fTG',
                  subplot(4,4,find(trials==tr));
                case 'mR'
                  subplot(4,4,find(trials==tr)+12);
                case 'lTG'
                  subplot(4,4,find(trials==tr)+7);
            end
          end
          if groupSize==25 && shelterSpace==25
             switch groupProtocol
                case 'fTG',
                  subplot(3,5,find(trials==tr));
                case 'mR'
                  subplot(3,5,find(trials==tr)+3);
                case 'lTG'
                  subplot(3,5,find(trials==tr)+10);
            end  
          end
        
        
        
        
        
        
