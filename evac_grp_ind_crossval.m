%% CROSS-VALIDATION 

%% set training and holdout trials; set shelter space

ss = 25;
[~,~,~,~,missing] = load_evac_data(0);
train_trials = trial_ix('ind',ss,1,missing);

makeFigs = true;  % allows any figures
legend1 = false;  % keep this off...
sumFigs = false;  % overrides full plots to make plots with ONE hit/miss
makeFigs_fit = false; % plots individual best fit model strategy
toSave = false;

% for cross-validation on group trials:
gPs = {'fTG','lTG','mR'};
gSs = [5,25];

% set parameters for sims
p = 100;
shelterSpace = ss;

% loop over group sizes and protocols
for groupSize = gSs
    
  nextNo = next_fig;
  for groupProtocol = gPs
    groupProtocol = char(groupProtocol);
    
    % set test trials and protocol strings
    test_trials = trial_ix(groupProtocol,ss,groupSize,missing);
    switch groupProtocol
                case 'fTG',
                  gpstr = 'FTG';
                case 'mR'
                  gpstr = 'MV';
                case 'lTG'
                  gpstr = 'LTG';
    end

    % load necessary data
    bins = 0:0.1:1.1;   % this rounds DOWN
    [Htrain,Jtrain,Theta,AvgCumEvac,missing,Q1,T1,P1,C1] = load_evac_data(0,train_trials,bins);
    [Htest,Jtest] = load_evac_data(0,test_trials,bins);
    [ntrials,nts] = size(Q1);
    [N,~] = size(P1);
    
    trEnd = zeros(1,size(Q1,1)); % we will store the time at which the trials ended
    for j = 1:size(Q1,1) % iterate through each trial
        trEnd(j) = find(Q1(j,:)==Q1(j,end),1,'first');
    end

    % first, perform (individual) fit on training trials with MLE

    % qform = @(Phit_,pv_) pv_(1).*(Phit_.^pv_(3))./(pv_(2).^pv_(3)+Phit_.^pv_(3));
    % Phit = 0:0.1:1;
    % qfit = @(pv_) qform(Phit(:,2:end),pv_);
    % MLfit = @(pvec) -1*sum(((sum(H(:,2:end))-sum(J(:,2:end))).*log(1-qfit(pvec)) + sum(J(:,2:end)).*log(qfit(pvec))));
    % startp = [1;0.5;10];

    qform = @(Phit_,pv_) pv_(1).*Phit_.^pv_(2);
    Phit = 0:0.1:1;
    qfit = @(pv_) qform(Phit(:,2:end-1),pv_);
    MLfit = @(pvec_) -1*sum(((sum(Htrain(:,2:end-1))-sum(Jtrain(:,2:end-1))).*log(1-qfit(pvec_)) + sum(Jtrain(:,2:end-1)).*log(qfit(pvec_))));
    startp = [1.1;1];
    options = optimoptions(@fminunc,'MaxFunEvals',10000);
    [params1,MLval] = fminunc(MLfit,startp,options);
    A = zeros(length(startp));
    b = zeros(length(startp),1);
    % constrains all prob. as a function of Phit to be <=1
    q_con = @(pvec_)(qform(Phit,pvec_)-1); 
    [theta_con,MLval_con] = fmincon(MLfit,startp,A,b,A,b,0,100,@(pvec_)q_add_eq(pvec_,q_con));

    if MLval<MLval_con && all(qform(0:0.1:1,params1)<=1)
        params = params1;
    else
        params = theta_con;
    end
    
% MAKE bestfit FIGURE if makeFigs = true
    if makeFigs_fit
        figure(1); hold all; plot(Phit,qform(Phit,params));
        xlabel('Phit'); ylabel('evacuation probability'); title('best-fit hill model, ind. trials, shelterSpace=50');
        %legend(['\Lambda = ' num2str(params(1)) ', k = ' num2str(params(2)) ', n = ' num2str(params(3))]);
        legend(['\alpha = ' num2str(params(1)) ', \beta = ' num2str(params(2))]);
        hold off;
    end

% Given training params, perform test on holdout trials
    Phit = 0:0.1:1;
    q = @(Phit_) qform(Phit_,params);              % model with params to test
    tot_err = zeros(numel(test_trials),1);
    ts_err = NaN*ones(numel(test_trials),nts);
    
    % preallocate for sims
    nevac_b = zeros(numel(test_trials),1);
    nevac_t = NaN*ones(numel(test_trials),nts);
    grp_dtime = zeros(numel(test_trials),N/groupSize,p);
    grp_grpIDs = zeros(numel(test_trials),p,N);
    grp_dtime_var = zeros(numel(test_trials),N/groupSize,p);
    grp_dtime_dist_samp = zeros(numel(test_trials),N/groupSize,p);
    grp_dtime_dist_empi = zeros(numel(test_trials),N/groupSize,p);

    %loop over test trials
    for tr=test_trials';
    disp(['test trial ' num2str(tr)]);
    [tr_hit,tr_miss] = tr_rep(ss,groupSize,groupProtocol);
    %if sumFigs && ~any([tr_hit,tr_miss]==tr)
    %    continue
    %end
    [T, P] = solve_master_binom_ss(N,q,Q1(tr,:),Phit,'ss',ss);
    Cexp = sum(bsxfun(@times,P,0:1:N),2);
    tbins = 0:1:nts-1;
    Cix = zeros(numel(tbins),1);
    for tx=tbins
      try  
        Cix(tbins==tx) = find(T>=tx,1,'first');
      catch
        Cix(tbins==tx) = Cix(find(tbins==tx)-1);
      end
    end
    tot_err(test_trials==tr) = Cexp(trEnd(tr)) - C1(tr,trEnd(tr));
    ts_err(test_trials==tr,1:trEnd(tr)) = Cexp(Cix(1:trEnd(tr)))' - C1(tr,1:trEnd(tr));

    err_pos = (sum(cumsum(P')<0.95));
    err_neg = (sum(cumsum(P')<0.05));
    stdevs = [err_pos-Cexp';Cexp'-err_neg];

    % randomly generate N individual evac times, p simulations
    % for comments see similar section above
    q = @(Phit_) qform(Phit_,params);
    tbins_s = 0:1:60;
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
    Tlen = endx + 1 - (endx==length(Q1(tr,:)));  % duration of simulation
    [ix,jx] = ind2sub(size(aix),find(~~aix));
    sx = sub2ind(size(aix),find(~~aix));
    nu = sparse(jx,1:length(sx),-1,N+1,length(sx)) + sparse(ix,1:length(sx),1);
    
    % this generates the samples with Gillespie algorithm...
    [t_mat,n_mat,r_mat,ns_mat,br_mat] = gillespie_tdep(theta,a,nu,P0,t0,Tlen,tbins_s,p);
    
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
    [Cgrp,tgrp,Cbin] = cum_evac(t_evac,grpIDs,groupProtocol,tbins_s,'ss',shelterSpace);
    
    % for each relevant trial, compute nevac_b and nevac(t)
    nevac_t(test_trials==tr,1:trEnd(tr)) = mean(Cbin(2:trEnd(tr)+1,:),2)' - C1(tr,1:trEnd(tr));
    nevac_b(test_trials==tr) = mean(Cbin(trEnd(tr)+1,:)) - C1(tr,trEnd(tr));
    
% MAKE comparison FIGURE if makeFigs = true    
if makeFigs    
    if legend1
        figure; %bcolor(P',T,0:1:N); colorbar; hold all; % probability distribution
        plot(T,Cexp,'k','LineWidth',2); hold all;  % expected value
        plot(tbins,C1(tr,:),'color',[0 0.5 0],'LineWidth',2);                         % empirical value
        plot(tbins,N*(floor((Q1(tr,:))*10)/10),'--','LineWidth',1);                 % Phit value
        confinterval=shadedErrorBar_2(T,Cexp,stdevs,{'color','k','LineWidth',2},1);
        set(confinterval.edge(1),'visible','off');
        set(confinterval.edge(2),'visible','off');
        xlabel('time'); ylabel('# evacuated'); title(['trial ' num2str(tr)]);
        axis([0 60 0 50]);
        legend('expected individual behavior','empirical group behavior','disaster trajectory','95% confidence');
        hold off;
    end

if sumFigs && ~any([tr_hit,tr_miss]==tr)
    continue
end    
figure(nextNo); hold all;
    if ~sumFigs
          if groupSize==5 && ss==50
            makeLeg = false;
            switch groupProtocol
                case 'fTG',
                  subplot(3,6,find(test_trials==tr));
                case 'mR'
                  subplot(3,6,[11,12,17,18]);
                  makeLeg = true;
                case 'lTG'
                  if find(test_trials==tr)<2  
                    subplot(3,6,find(test_trials==tr)+9);
                  else
                    subplot(3,6,find(test_trials==tr)+11);  
                  end
            end
          end
          if groupSize==25 && ss==50
             makeLeg = false;
             switch groupProtocol
                case 'fTG',
                  subplot(4,5,find(test_trials==tr));
                case 'mR'
                  if find(test_trials==tr)==1  
                  subplot(4,5,5);
                  else
                      if find(test_trials==tr)<7
                          subplot(4,5,find(test_trials==tr)+4);
                      else
                          subplot(4,5,[14,15,19,20]);
                          makeLeg = true;
                      end
                  end
                case 'lTG'
                  if find(test_trials==tr)<4
                    subplot(4,5,find(test_trials==tr)+10);
                  else
                    subplot(4,5,find(test_trials==tr)+12);
                  end
            end  
          end
          if groupSize==5 && ss==25
            makeLeg = false;
            switch groupProtocol
                case 'fTG',
                  subplot(3,7,find(test_trials==tr));
                case 'mR'
                  if find(test_trials==tr)<6
                    subplot(3,7,find(test_trials==tr)+7);
                  else
                    makeLeg = true;
                    subplot(3,7,[13,14,20,21]);
                  end
                case 'lTG'
                  subplot(3,7,find(test_trials==tr)+14);
            end
          end
          if groupSize==25 && ss==25
             makeLeg = false;
             switch groupProtocol
                case 'fTG',
                  subplot(4,5,find(test_trials==tr));
                case 'mR'
                  if find(test_trials==tr)<3  
                    subplot(4,5,find(test_trials==tr)+3);
                  else
                    subplot(4,5,find(test_trials==tr)+3);  
                  end
                case 'lTG'
                  if find(test_trials==tr)<4  
                    subplot(4,5,find(test_trials==tr)+10);
                  else
                    if find(test_trials==tr)<7  
                      subplot(4,5,find(test_trials==tr)+12);
                    else
                      makeLeg = true;
                      subplot(4,5,[14,15,19,20]);  
                    end
                  end
            end  
          end
    end
    if sumFigs && any([tr_hit,tr_miss]==tr)  % for sumFigs
         horm = Q1(tr,end);  % (=0 if miss, =1 if hit)
         assert(find([tr_hit,tr_miss]==tr)==(2-horm));
         makeLeg = false;
         subplot(3,2,2*find(ismemvar(gPs,groupProtocol))-horm);
         if horm && strcmp(groupProtocol,'mR')
            makeLeg = true;
         end
    end  % end fig size/arrangement preamble
          
%figure; %bcolor(P',T,0:1:N); colorbar; hold all; % probability distribution
        plot(tbins,N*(floor((Q1(tr,:))*10)/10),'--','color',[0 0.5 0],'LineWidth',1);  hold all;
        plot(tbins,C1(tr,:),'k','LineWidth',2.5);                         % empirical value
        plot(T,Cexp,'color',[0.5 0 0.5],'LineWidth',2); hold all;  % expected value
        plot(tbins_s,mean(Cbin,2),'b','LineWidth',2); % simulated value
        confinterval=shadedErrorBar_2(T,Cexp,stdevs,{'color',[0.5 0 0.5],'LineWidth',2},1);
        set(confinterval.edge(1),'visible','off');
        set(confinterval.edge(2),'visible','off');
        for i=1:10%p
                plot([0;tgrp(:,i)],[0;Cgrp(:,i)],':','LineWidth',1); 
                %plot(t_mat(r_mat(:,i)~=0,i),ix(r_mat(r_mat(:,i)~=0,i))-1,':o'); hold all;
        end
        xlabel('time'); ylabel('# evacuated'); title(['trial ' num2str(tr) ', ' gpstr]);
        axis([0 60 0 50]);
        if makeLeg
            legend('disaster trajectory','empirical group behavior','expected naive cross-val behavior',...%'95% confidence',
                   'mean simulated behavior');%,'example simulation instances');
        end
        %suptitle(['Shelter Capacity ' num2str(ss) ', Groups of ' num2str(groupSize)]);
        title(trial_conv(tr,ss,groupSize,groupProtocol));
        hold off; 
end  % end makeFigs

    end  % end loop over test trials

% SAVE the relevant errors for crossval AND sims
savingstr = [lower(gpstr),num2str(groupSize)];
eval(['rmse_' savingstr ' = rmse(ts_err'');']);
eval(['tse_' savingstr ' = ts_err;']);
eval(['rmse_' savingstr 's = rmse(nevac_t'');']);
eval(['tse_' savingstr 's = nevac_t;']);

  end  % end loop over protocol
  
  %close all;  % clear current figs
  
end  % end loop over group size

if toSave
save(['ind_grp_errs_ss' num2str(ss) 't.mat'],'rmse_ftg5',...
                                                 'rmse_ftg25',...
                                                 'rmse_ltg5',...
                                                 'rmse_ltg25',...
                                                 'rmse_mv5',...
                                                 'rmse_mv25',...
                                                 'tse_ftg5',...
                                                 'tse_ftg25',...
                                                 'tse_ltg5',...
                                                 'tse_ltg25',...
                                                 'tse_mv5',...
                                                 'tse_mv25',...
                                                 'rmse_ftg5s',...
                                                 'rmse_ftg25s',...
                                                 'rmse_ltg5s',...
                                                 'rmse_ltg25s',...
                                                 'rmse_mv5s',...
                                                 'rmse_mv25s',...
                                                 'tse_ftg5s',...
                                                 'tse_ftg25s',...
                                                 'tse_ltg5s',...
                                                 'tse_ltg25s',...
                                                 'tse_mv5s',...
                                                 'tse_mv25s',...
                                                 '-v7.3');
end
%% make RMSE plot for a particular ss and groupSize

ss = 25;
nextNo = next_fig;
add_ind = true;
add_ind_1 = false;

%try_labels = {sprintf('FTG LTG MV\nCV'),sprintf('FTG LTG MV\nSIM')};
for groupSize = gSs
prots = {'ftg','ltg','mv'};
load(['ind_grp_errs_ss' num2str(ss) 't.mat']);
if add_ind
    load(['ind_errs_ss' num2str(ss) 't.mat']);
end
                  
% separate by hit/miss rmse_dat
rmse_dat = zeros(2,3,2); % (cv/sim, protocols, hit/miss)
for gg = gPs
    gg = char(gg);
    tr_horm = Q1(trial_ix(gg,ss,groupSize,missing),end);
    ptx = find(strcmp(gPs,gg));
    if numel(tr_horm)
      % hits CV
      eval(char(strcat('rmse_dat(1,ptx,2) = mean(rmse_',prots(ptx),num2str(groupSize),'(~~tr_horm))*sum(~~tr_horm(:))/numel(tr_horm);')));
      % misses CV
      eval(char(strcat('rmse_dat(1,ptx,1) = mean(rmse_',prots(ptx),num2str(groupSize),'(~tr_horm))*sum(~tr_horm(:))/numel(tr_horm);')));
      % hits SIM
      eval(char(strcat('rmse_dat(2,ptx,2) = mean(rmse_',prots(ptx),num2str(groupSize),'s(~~tr_horm))*sum(~~tr_horm(:))/numel(tr_horm);')));
      % misses SIM
      eval(char(strcat('rmse_dat(2,ptx,1) = mean(rmse_',prots(ptx),num2str(groupSize),'s(~tr_horm))*sum(~tr_horm(:))/numel(tr_horm);')));
    end
end

if add_ind
    rmse_ind_dat = NaN*ones(1,3,2);
    tr_horm = Q1(trial_ix('ind',ss,1,missing),end);
    if numel(tr_horm)
        rmse_ind_dat(1,2,2) = mean(rmse_ind(~~tr_horm))*sum(~~tr_horm(:))/numel(tr_horm);  % hits CV
        rmse_ind_dat(1,2,1) = mean(rmse_ind(~tr_horm))*sum(~tr_horm(:))/numel(tr_horm);  % misses CV
    end
    rmse_dat = [rmse_ind_dat;rmse_dat];
    try_labels = {'Ind','Grp CV','Grp SIM'};
end

if add_ind_1
    rmse_ind_dat = zeros(2,1,2);
    tr_horm = Q1(trial_ix('ind',ss,1,missing),end);
    if numel(tr_horm)
        rmse_ind_dat(1,1,2) = mean(rmse_ind(~~tr_horm))*sum(~~tr_horm(:))/numel(tr_horm);  % hits CV
        rmse_ind_dat(1,1,1) = mean(rmse_ind(~tr_horm))*sum(~tr_horm(:))/numel(tr_horm);  % misses CV
        rmse_ind_dat(2,1,2) = NaN;  % hits sim
        rmse_ind_dat(2,1,1) = NaN;  % misses sim
    end
    rmse_dat = [rmse_ind_dat,rmse_dat];
end

figure(nextNo); hold on;
        subplot(1,2,find(gSs==groupSize));
        plotBarStackGroups(rmse_dat,try_labels); 
        legend('miss trials','hit trials');
        title(['Shelter Space = ' num2str(ss) ', Group Size = ' num2str(groupSize)]);
        ylabel('RMS error (evacuations), mean over trials');
        hold off;
end

%% summary figs....

% for a particular shelter capacity
ss = 25;

gSs = [5,25];
tbins = 0:1:nts-1;
for groupSize = gSs;
    nextNo = next_fig;
    gix = find(gSs==groupSize);

        strg1 = 'Naive Cross-validation Error';
        strg2 = 'Grouped Individual Simulation Error';
        load(['ind_grp_crossval_ss' num2str(ss) 't.mat']); % t suffix = crossval; s suffix = sims
  figure;
    subplot(2,1,1);
    eval(['plot(tbins,mean(tse_ftg' num2str(groupSize) ')'',''b'',''LineWidth'',2);']); hold all;
    eval(['plot(tbins,mean(tse_ltg' num2str(groupSize) ')'',''m'',''LineWidth'',2);']); 
    eval(['tsmr=tse_mv' num2str(groupSize) ';']);
    if numel(tsmr)==length(tsmr)
      eval(['plot(tbins,tse_mv' num2str(groupSize) ''',''g'',''LineWidth'',2);']);
    else
      eval(['plot(tbins,mean(tse_mv' num2str(groupSize) ')'',''g'',''LineWidth'',2);']);
    end
    eval(['plot(tbins,tse_ftg' num2str(groupSize) ''',''b:'');']); hold all;
    eval(['plot(tbins,tse_ltg' num2str(groupSize) ''',''m:'');']); 
    eval(['plot(tbins,tse_mv' num2str(groupSize) ''',''g:'');']);
    xlabel('time step'); ylabel('prediction error in total evacuations');
    legend('FTG','LTG','MV'); axis([0 nts -ss ss]);
    %title(['shelter capacity ' num2str(ss) ', group size ' num2str(groupSize)]);
    title(strg1);
    
    subplot(2,1,2);
    eval(['plot(tbins,nanmean(tse_ftg' num2str(groupSize) 's)'',''b'',''LineWidth'',2);']); hold all;
    eval(['plot(tbins,nanmean(tse_ltg' num2str(groupSize) 's)'',''m'',''LineWidth'',2);']); 
    eval(['tsmr=tse_mv' num2str(groupSize) ';']);
    if numel(tsmr)==length(tsmr)
      eval(['plot(tbins,tse_mv' num2str(groupSize) 's'',''g'',''LineWidth'',2);']);
    else
      eval(['plot(tbins,nanmean(tse_mv' num2str(groupSize) 's)'',''g'',''LineWidth'',2);']);
    end
    eval(['plot(tbins,tse_ftg' num2str(groupSize) 's'',''b:'');']); hold all;
    eval(['plot(tbins,tse_ltg' num2str(groupSize) 's'',''m:'');']); 
    eval(['plot(tbins,tse_mv' num2str(groupSize) 's'',''g:'');']);
    xlabel('time step'); ylabel('prediction error in total evacuations');
    legend('FTG','LTG','MV'); axis([0 nts -ss ss]);
    %title(['shelter capacity ' num2str(ss) ', group size ' num2str(groupSize)]);
    title(strg2);
end

%% RUN LOOCV

ss = 25;
[~,~,~,~,missing] = load_evac_data(0);
ind_trials = trial_ix('ind',ss,1,missing);

rmse_ind = zeros(numel(ind_trials),1);
tse_ind = zeros(numel(ind_trials),nts);

% for LOOCV on ind trials:
for test_trials = ind_trials(:)'
    train_trials = removeval(ind_trials,test_trials);
    
    % load necessary data
    bins = 0:0.1:1.1;   % this rounds DOWN
    [Htrain,Jtrain,Theta,AvgCumEvac,missing,Q1,T1,P1,C1] = load_evac_data(0,train_trials,bins);
    [Htest,Jtest] = load_evac_data(0,test_trials,bins);
    [ntrials,nts] = size(Q1);
    [N,~] = size(P1);
    plotH = test_trials(find(Q1(test_trials,end)==1,1,'first'));
    plotM = test_trials(find(Q1(test_trials,end)==0,1,'first'));

    trEnd = zeros(1,size(Q1,1)); % we will store the time at which the trials ended
    for j = 1:size(Q1,1) % iterate through each trial
        trEnd(j) = find(Q1(j,:)==Q1(j,end),1,'first');
    end

    % first, perform (individual) fit on training trials with MLE

    % qform = @(Phit_,pv_) pv_(1).*(Phit_.^pv_(3))./(pv_(2).^pv_(3)+Phit_.^pv_(3));
    % Phit = 0:0.1:1;
    % qfit = @(pv_) qform(Phit(:,2:end),pv_);
    % MLfit = @(pvec) -1*sum(((sum(H(:,2:end))-sum(J(:,2:end))).*log(1-qfit(pvec)) + sum(J(:,2:end)).*log(qfit(pvec))));
    % startp = [1;0.5;10];

    qform = @(Phit_,pv_) pv_(1).*Phit_.^pv_(2);
    Phit = 0:0.1:1;
    qfit = @(pv_) qform(Phit(:,2:end-1),pv_);
    MLfit = @(pvec_) -1*sum(((sum(Htrain(:,2:end-1))-sum(Jtrain(:,2:end-1))).*log(1-qfit(pvec_)) + sum(Jtrain(:,2:end-1)).*log(qfit(pvec_))));
    startp = [1.1;1];
    options = optimoptions(@fminunc,'MaxFunEvals',10000);
    [params1,MLval] = fminunc(MLfit,startp,options);
    A = zeros(length(startp));
    b = zeros(length(startp),1);
    % constrains all prob. as a function of Phit to be <=1
    q_con = @(pvec_)(qform(Phit,pvec_)-1); 
    [theta_con,MLval_con] = fmincon(MLfit,startp,A,b,A,b,0,100,@(pvec_)q_add_eq(pvec_,q_con));

    if MLval<MLval_con && all(qform(0:0.1:1,params1)<=1)
        params = params1;
    else
        params = theta_con;
    end

    % Given training params, perform test on holdout trials
    Phit = 0:0.1:1;
    q = @(Phit_) qform(Phit_,params);              % model with params to test
    tot_err = zeros(numel(test_trials),1);
    ts_err = NaN*ones(numel(test_trials),nts);
    
    %loop over test trials
    for tr=test_trials';
    disp(['test trial ' num2str(tr)]);
    
    [T, P] = solve_master_binom_ss(N,q,Q1(tr,:),Phit,'ss',ss);
    Cexp = sum(bsxfun(@times,P,0:1:N),2);
    tbins = 0:1:nts-1;
    Cix = zeros(numel(tbins),1);
    for tx=tbins
      try  
        Cix(tbins==tx) = find(T>=tx,1,'first');
      catch
        Cix(tbins==tx) = Cix(find(tbins==tx)-1);
      end
    end
    tot_err(test_trials==tr) = Cexp(trEnd(tr)) - C1(tr,trEnd(tr));
    ts_err(test_trials==tr,1:trEnd(tr)) = Cexp(Cix(1:trEnd(tr)))' - C1(tr,1:trEnd(tr));

    err_pos = (sum(cumsum(P')<0.95));
    err_neg = (sum(cumsum(P')<0.05));
    stdevs = [err_pos-Cexp';Cexp'-err_neg];
    
    end
    
    rmse_ind(ind_trials==test_trials) = rmse(ts_err');
    tse_ind(ind_trials==test_trials,:) = ts_err';
    
end

	% SAVE the relevant errors for crossval AND sims
    save(['ind_errs_ss' num2str(ss) 't.mat'],'rmse_ind',...
                                             'tse_ind',...
                                             '-v7.3');
    


%% MAKE SUMMARY PLOTS (OLD)


% for a particular shelter capacity
ss_list = 25;
for ss = ss_list;
figure(find(ss_list==ss));
gSs = [5,25];
gSs = gSs(gSs<=ss);
tbins = 0:1:nts-1;
for groupSize = gSs;
    gix = find(gSs==groupSize);
  for tp=[1,2];  
    if tp==1
        strg = 'Naive Cross-validation Error';
        load(['ind_grp_crossval_ss' num2str(ss) 't.mat']); % t suffix = crossval; s suffix = sims
    else
        strg = 'Grouped Individual Simulation Error';
        load(['ind_grp_crossval_ss' num2str(ss) 's.mat']); % t suffix = crossval; s suffix = sims
    end
    subplot(2,2,2*gix-(2-tp));
    eval(['plot(tbins,mean(tse_ftg' num2str(groupSize) ')'',''b'',''LineWidth'',2);']); hold all;
    eval(['plot(tbins,mean(tse_ltg' num2str(groupSize) ')'',''m'',''LineWidth'',2);']); 
    eval(['tsmr=tse_mr' num2str(groupSize) ';']);
    if numel(tsmr)==length(tsmr)
      eval(['plot(tbins,tse_mr' num2str(groupSize) ''',''g'',''LineWidth'',2);']);
    else
      eval(['plot(tbins,mean(tse_mr' num2str(groupSize) ')'',''g'',''LineWidth'',2);']);
    end
    eval(['plot(tbins,tse_ftg' num2str(groupSize) ''',''b--'');']); hold all;
    eval(['plot(tbins,tse_ltg' num2str(groupSize) ''',''m--'');']); 
    eval(['plot(tbins,tse_mr' num2str(groupSize) ''',''g--'');']);
    xlabel('time step'); ylabel('prediction error in total evacuations');
    legend('FTG','LTG','MV'); axis([0 nts -ss ss]);
    %title(['shelter capacity ' num2str(ss) ', group size ' num2str(groupSize)]);
    title([strg ', group size ' num2str(groupSize)]);
  end
end
end

%% old
% for a particular shelter capacity
ss_list = [50,25,5];
figure(3);
for ss = ss_list;
gSs = [5,25];
gSs = gSs(gSs<=ss);
load(['ind_grp_crossval_ss' num2str(ss) 't.mat']); % t suffix = crossval; s suffix = sims
load(['ind_grp_crossval_ss' num2str(ss) 's.mat']); % t suffix = crossval; s suffix = sims
for groupSize = gSs;
    gix = find(gSs==groupSize);
    %subplot(1,numel(sizes),gix);
    subplot(3,2,find(ss_list==ss)*2+gix-2);
    eval(['plot(tbins,mean(tse_ftg' num2str(groupSize) ')'',''b'',''LineWidth'',2);']); hold all;
    eval(['plot(tbins,mean(tse_ltg' num2str(groupSize) ')'',''m'',''LineWidth'',2);']); 
    eval(['plot(tbins,mean(tse_mr' num2str(groupSize) ')'',''g'',''LineWidth'',2);']);
    eval(['plot(tbins,tse_ftg' num2str(groupSize) ''',''b--'');']); hold all;
    eval(['plot(tbins,tse_ltg' num2str(groupSize) ''',''m--'');']); 
    eval(['plot(tbins,tse_mr' num2str(groupSize) ''',''g--'');']);
    xlabel('time step'); ylabel('RMS prediction error in cum. evacuations');
    legend('FTG','LTG','MV'); axis([0 nts -ss ss]);
    title(['shelter capacity ' num2str(ss) ', group size ' num2str(groupSize)]);
end
end

%% MAKE PLOTS (old)

% for a particular shelter capacity

ss = 50;
gSs = [5,25];
gSs = gSs(gSs<=ss);
load(['ind_grp_crossval_ss' num2str(ss) '.mat']);
figure; 
for groupSize = gSs;
    gix = find(gSs==groupSize);
    subplot(1,numel(gSs),gix);
    eval(['plot(tbins,RMSE_ftg' num2str(groupSize) ');']); hold all;
    eval(['plot(tbins,RMSE_ltg' num2str(groupSize) ');']); 
    eval(['plot(tbins,RMSE_mr' num2str(groupSize) ');']); 
    xlabel('time step'); ylabel('RMS prediction error in cum. evacuations');
    legend('FTG','LTG','MV'); axis([0 nts 0 ss]);
    title(['shelter capacity ' num2str(ss) ', group size ' num2str(groupSize)]);
end
