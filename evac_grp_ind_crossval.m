%% Comparing Group to Individual behavior w CROSS-VALIDATION and SIMULATION
% this code generates parts A, B, and C of group-individual comparison
% figures in decision modelling paper
%cd kimberly/Documents/evacuation
addpath(genpath(pwd));

%% set training and holdout trials; set shelter space

ss = 25;
[~,~,~,~,missing] = load_evac_data(0);
train_trials = trial_ix('ind',ss,1,missing);

makeFigs = true;  % allows any figures
sumFigs = true;  % overrides full plots to make plots with ONE example trial for each hit/miss
                  % NOTE if sumFigs = false && makeFigs = true, plots
                  % will be made with ALL trials included
makeFigs_fit = false; % plots individual best fit model strategy
toSave = false;   % if toSave = true, ALL simulations will be run and all errors saved

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
    [Htrain,Jtrain,Theta,AvgCumEvac,missing,Q1,~,P1,C1] = load_evac_data(0,train_trials,bins);
    [Htest,Jtest] = load_evac_data(0,test_trials,bins);
    [ntrials,n_ts] = size(Q1);
    [N,~] = size(P1);
    
    trEnd = zeros(1,size(Q1,1)); % we will store the time at which the trials ended
    for j = 1:size(Q1,1) % iterate through each trial
        trEnd(j) = find(Q1(j,:)==Q1(j,end),1,'first');
    end

    % first, perform (individual) fit on training trials with maximum likelihood estimation
    qform = @(Phit_,pv_) pv_(1).*Phit_.^pv_(2);
    %startp = [1.1;1];
    startp = [0.9;10];
    Phit = 0:0.1:1;
    params = ML_fit_beta(qform,Phit(:,2:end),Htrain(:,2:end),Jtrain(:,2:end),startp);

    % MAKE bestfit FIGURE if makeFigs_fit = true
    if makeFigs_fit
        figure(next_fig+1); hold all; plot(Phit,qform(Phit,params));
        xlabel('Phit'); ylabel('evacuation probability'); title('best-fit hill model, ind. trials, shelterSpace=50');
        %legend(['\Lambda = ' num2str(params(1)) ', k = ' num2str(params(2)) ', n = ' num2str(params(3))]);
        legend(['\alpha = ' num2str(params(1)) ', \beta = ' num2str(params(2))]);
        hold off;
    end

    % preallocate for cross-validation error
    tot_err = zeros(numel(test_trials),1);
    ts_err = NaN*ones(numel(test_trials),n_ts);
    
    % preallocate for simulation error
    nevac_t = NaN*ones(numel(test_trials),n_ts);
    
    % model with parameters to test
	q = @(Phit_) qform(Phit_,params); 
    
    % loop over test trials
    for tr=test_trials';
    disp(['test trial ' num2str(tr)]);
    
        % find the example trials to display
        [tr_hit,tr_miss] = tr_rep(ss,groupSize,groupProtocol);
        
        % if ONLY sumFigs, skip calculations for non-example trials
        if sumFigs && ~toSave && ~any([tr_hit,tr_miss]==tr)
            continue;
        end

        % solve the master equation with best-fit strategy
        [T, P] = solve_master_binom_ss(N,q,Q1(tr,:),Phit,'ss',ss);
        Cexp = sum(bsxfun(@times,P,0:1:N),2);  % expected value

        % bin the expected value of master equation solution into time steps
        tbins = 0:1:n_ts-1;
        Cix = zeros(numel(tbins),1);
        for tx=tbins
          try  
            Cix(tbins==tx) = find(T>=tx,1,'first');
          catch
            Cix(tbins==tx) = Cix(find(tbins==tx)-1);
          end
        end

        % record total cross-validation error and error on each time step
        tot_err(test_trials==tr) = Cexp(trEnd(tr)) - C1(tr,trEnd(tr));
        ts_err(test_trials==tr,1:trEnd(tr)) = Cexp(Cix(1:trEnd(tr)))' - C1(tr,1:trEnd(tr));

        % find 3-sigma envelopes
        err_pos = (sum(cumsum(P')<0.997));
        err_neg = (sum(cumsum(P')<0.003));
        stdevs = [err_pos-Cexp';Cexp'-err_neg];

        % Simulation: prepare inputs
        q = @(Phit_) qform(Phit_,params);  % q as a function of Phit for simulation
        tbins_sim = 0:1:60;
        a = @(t_,P_,theta_) a_vec(t_,P_,theta_{1},theta_{2}); % this returns a propensity vector of only lower-triangular vals

        theta = {makeA(N,q(Phit),ss); Q1(tr,:)}; % these are the non-explicitly-time-dependent input parameters

        P0 = [1;zeros(N,1)];   % this is the initial condition
        t0 = 0;                % starting time of simulation
        endx = trEnd(tr);
        Tlen = endx + 1 - (endx==length(Q1(tr,:)));  % duration of simulation

        aix = tril(ones(N+1),-1);  % locations of lower-triangular entries
        [ix,jx] = ind2sub(size(aix),find(~~aix)); % x and y indices
        sx = sub2ind(size(aix),find(~~aix)); % overall indices
        nu = sparse(jx,1:length(sx),-1,N+1,length(sx)) + sparse(ix,1:length(sx),1); % transition instructions for all possible transitions

        % this generates the samples with Gillespie algorithm...
        % p total simulations, N evacuation times for each
        [t_mat,n_mat,r_mat,ns_mat,br_mat] = gillespie_tdep(theta,a,nu,P0,t0,Tlen,tbins_sim,p);

        % this creates a list of evac times & group IDs for each participant
        [t_evac,grpIDs] = group_evac_times(N,p,groupSize,t_mat,r_mat,ix,jx);

        % this creates cumulative evac plots given group IDs and evac times
        [Cgrp,tgrp,Cbin] = cum_evac(t_evac,grpIDs,groupProtocol,tbins_sim,'ss',shelterSpace);

        % for each relevant trial, compute nevac_t (# evacuated error at each sim timestep)
        nevac_t(test_trials==tr,1:trEnd(tr)) = mean(Cbin(2:trEnd(tr)+1,:),2)' - C1(tr,1:trEnd(tr));
    
     if makeFigs    
      
      % if summary figures are called for and this is not an example trial,
      % do not make a figure here
      if sumFigs && ~any([tr_hit,tr_miss]==tr)
        continue
      end
      
      % else, create the figure with the number chosen for this groupSize
      figure(nextNo); hold all;
        % set figure format & subplot location
        if ~sumFigs  % for no summary figures, use fig format with many panels
            makeLeg = fig_format_many(tr,test_trials,ss,groupProtocol,groupSize);
        end
        if sumFigs && any([tr_hit,tr_miss]==tr)  % for sumFigs, use summary format
            makeLeg = fig_format_sum(tr,Q1(tr,end),tr_hit,tr_miss,gPs,groupProtocol);
        end  

        % make plot
            plot(tbins,N*(floor((Q1(tr,:))*10)/10),'--','LineWidth',1); hold all;% Phit value
            plot(tbins,C1(tr,:),'color','k','LineWidth',2);  % empirical value
            plot(T,Cexp,'color',[0.5 0 0.5],'LineWidth',2); % expected value
            plot(tbins_sim,mean(Cbin,2),'color',[0 0.5 0],'LineWidth',2);
            confinterval=shadedErrorBar_2(T,Cexp,stdevs,{'color',[0.5 0 0.5],'LineWidth',2},1);
            set(confinterval.edge(1),'visible','off');
            set(confinterval.edge(2),'visible','off');
            for i=1:10%p
                    plot([0;tgrp(:,i)],[0;Cgrp(:,i)],':','LineWidth',1); 
                    %plot(t_mat(r_mat(:,i)~=0,i),ix(r_mat(r_mat(:,i)~=0,i))-1,':o'); hold all;
            end
            xlabel('time'); ylabel('# evacuated'); title(['trial ' num2str(tr) ', ' gpstr]);
            axis([0 n_ts 0 N]);
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
gSs = [5,25];
gPs = {'fTG','lTG','mR'};
nextNo = next_fig;
add_ind = true;
add_ind_1 = false;
[~,~,~,~,missing,Q1,T1,P1,C1] = load_evac_data(0);

%try_labels = {sprintf('FTG LTG MV\nCV'),sprintf('FTG LTG MV\nSIM')};
for groupSize = gSs
prots = {'ftg','ltg','mv'};
load(['ind_grp_errs_ss' num2str(ss) 't.mat']);
if add_ind  % includes an individual-game bar
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

if add_ind  % use this to add individual errors to plot
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

%% summary error figures....

% for a particular shelter capacity
ss = 25;

gSs = [5,25];
tbins = 0:1:n_ts-1;
for groupSize = gSs;
    nextNo = next_fig;
    gix = find(gSs==groupSize);

        strg1 = 'Naive Cross-validation Error';
        strg2 = 'Grouped Individual Simulation Error';
        load(['IndGrp_Results/ind_grp_crossval_ss' num2str(ss) 't.mat']); % t suffix = crossval; s suffix = sims
  figure(nextNo);
    subplot(2,1,1);
    eval(['plot(tbins,nanmean(tse_ftg' num2str(groupSize) ')'',''b'',''LineWidth'',2);']); hold all;
    eval(['plot(tbins,nanmean(tse_ltg' num2str(groupSize) ')'',''m'',''LineWidth'',2);']); 
    eval(['tsmr=tse_mv' num2str(groupSize) ';']);
    if numel(tsmr)==length(tsmr)
      eval(['plot(tbins,tse_mv' num2str(groupSize) ''',''g'',''LineWidth'',2);']);
    else
      eval(['plot(tbins,nanmean(tse_mv' num2str(groupSize) ')'',''g'',''LineWidth'',2);']);
    end
    eval(['plot(tbins,tse_ftg' num2str(groupSize) ''',''b:'');']); hold all;
    eval(['plot(tbins,tse_ltg' num2str(groupSize) ''',''m:'');']); 
    eval(['plot(tbins,tse_mv' num2str(groupSize) ''',''g:'');']);
    eval(['shdF = shadedErrorBar(tbins,nanmean(tse_ftg' num2str(groupSize) ',1),nanstd(tse_ftg' num2str(groupSize) ',0,1),{''color'',''b''},1,1);']);
    set(shdF.edge(1),'visible','off'); set(shdF.edge(2),'visible','off');
    eval(['shdL = shadedErrorBar(tbins,nanmean(tse_ltg' num2str(groupSize) ',1),nanstd(tse_ltg' num2str(groupSize) ',0,1),{''color'',''m''},1,1);']);
    set(shdL.edge(1),'visible','off'); set(shdL.edge(2),'visible','off');
    eval(['shdM = shadedErrorBar(tbins,nanmean(tse_mv' num2str(groupSize) ',1),nanstd(tse_mv' num2str(groupSize) ',0,1),{''color'',''g''},1,1);']);
    set(shdM.edge(1),'visible','off'); set(shdM.edge(2),'visible','off');
    xlabel('time step'); ylabel('prediction error in total evacuations');
    legend('FTG','LTG','MV'); axis([0 n_ts -ss ss]);
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
    eval(['shdFs = shadedErrorBar(tbins,nanmean(tse_ftg' num2str(groupSize) 's,1),nanstd(tse_ftg' num2str(groupSize) 's,0,1),{''color'',''b''},1,1);']);
    set(shdFs.edge(1),'visible','off'); set(shdFs.edge(2),'visible','off');
    eval(['shdLs = shadedErrorBar(tbins,nanmean(tse_ltg' num2str(groupSize) 's,1),nanstd(tse_ltg' num2str(groupSize) 's,0,1),{''color'',''m''},1,1);']);
    set(shdLs.edge(1),'visible','off'); set(shdLs.edge(2),'visible','off');
    eval(['shdMs = shadedErrorBar(tbins,nanmean(tse_mv' num2str(groupSize) 's,1),nanstd(tse_mv' num2str(groupSize) 's,0,1),{''color'',''g''},1,1);']);
    set(shdMs.edge(1),'visible','off'); set(shdMs.edge(2),'visible','off');
    xlabel('time step'); ylabel('prediction error in total evacuations');
    legend('FTG','LTG','MV'); axis([0 n_ts -ss ss]);
    %title(['shelter capacity ' num2str(ss) ', group size ' num2str(groupSize)]);
    title(strg2);
end

%% RUN LOOCV

ss = 25;
[~,~,~,~,missing] = load_evac_data(0);
ind_trials = trial_ix('ind',ss,1,missing);

rmse_ind = zeros(numel(ind_trials),1);
tse_ind = zeros(numel(ind_trials),n_ts);

% for LOOCV on ind trials:
for test_trials = ind_trials(:)'
    train_trials = removeval(ind_trials,test_trials);
    
    % load necessary data
    bins = 0:0.1:1.1;   % this rounds DOWN
    [Htrain,Jtrain,Theta,AvgCumEvac,missing,Q1,T1,P1,C1] = load_evac_data(0,train_trials,bins);
    [Htest,Jtest] = load_evac_data(0,test_trials,bins);
    [ntrials,n_ts] = size(Q1);
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
    ts_err = NaN*ones(numel(test_trials),n_ts);
    
    %loop over test trials
    for tr=test_trials';
    disp(['test trial ' num2str(tr)]);
    
    [T, P] = solve_master_binom_ss(N,q,Q1(tr,:),Phit,'ss',ss);
    Cexp = sum(bsxfun(@times,P,0:1:N),2);
    tbins = 0:1:n_ts-1;
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
    


