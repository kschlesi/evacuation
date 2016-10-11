%% CROSS-VALIDATION 

%% set training and holdout trials; set shelter space

ss = 25;
[~,~,~,~,missing] = load_evac_data(0);
train_trials = trial_ix('ind',ss,1,missing);

% for LOOCV on ind trials:
%loo_ix = 3;
%test_trials = train_trials(loo_ix);
%train_trials = removeval(train_trials,test_trials);

% for cross-validation on group trials:
groupProtocol = 'mR';
groupSize = 25;
test_trials = trial_ix(groupProtocol,ss,groupSize,missing);
switch groupProtocol
                case 'fTG',
                  gpstr = 'FTG';
                case 'mR'
                  gpstr = 'MV';
                case 'lTG'
                  gpstr = 'LTG';
end

makeFigs = true;
legend1 = false;

%% load necessary data

%bins = -0.05:0.1:1.05;  % this gives the bin EDGES for Phit values
bins = 0:0.1:1.1;   % this rounds DOWN
[Htrain,Jtrain,Theta,AvgCumEvac,missing,Q1,T1,P1,C1] = load_evac_data(0,train_trials,bins);
[Htest,Jtest] = load_evac_data(0,test_trials,bins);
[ntrials,nts] = size(Q1);
[N,~] = size(P1);

trEnd = zeros(1,size(Q1,1)); % we will store the time at which the trials ended
for j = 1:size(Q1,1) % iterate through each trial
    trEnd(j) = find(Q1(j,:)==Q1(j,end),1,'first');
end

%% first, perform (individual) fit on training trials
% with maximum likelihood estimation

% qform = @(Phit_,pv_) pv_(1).*(Phit_.^pv_(3))./(pv_(2).^pv_(3)+Phit_.^pv_(3));
% startp = [1;0.5;10];

qform = @(Phit_,pv_) pv_(1).*Phit_.^pv_(2);
startp = [1.1;1];
Phit = 0:0.1:1;
params = ML_fit_beta(qform,Phit(:,2:end-1),H(:,2:end-1),J(:,2:end-1),startp);

if makeFigs
figure(1); hold all; plot(Phit,qform(Phit,params));
xlabel('Phit'); ylabel('evacuation probability'); title('best-fit hill model, ind. trials, shelterSpace=50');
%legend(['\Lambda = ' num2str(params(1)) ', k = ' num2str(params(2)) ', n = ' num2str(params(3))]);
legend(['\alpha = ' num2str(params(1)) ', \beta = ' num2str(params(2))]);
hold off;
end

%% Given training params, perform test on holdout trials

Phit = 0:0.1:1;
q = @(Phit_) qform(Phit_,params);              % model with params to test
tot_err = zeros(numel(test_trials),1);
ts_err = zeros(numel(test_trials),nts);
%loop over test trials
for tr=test_trials'
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
        
figure(2); hold all;
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
          
    % make figure
        plot(T,Cexp,'k','LineWidth',2); hold all;  % expected value
        plot(tbins,C1(tr,:),'color',[0 0.5 0],'LineWidth',2);                         % empirical value
        plot(tbins,N*(floor((Q1(tr,:))*10)/10),'--','LineWidth',1);                 % Phit value
        confinterval=shadedErrorBar_2(T,Cexp,stdevs,{'color','k','LineWidth',2},1);
        set(confinterval.edge(1),'visible','off');
        set(confinterval.edge(2),'visible','off');
        xlabel('time'); ylabel('# evacuated'); title(['trial ' num2str(tr) ', ' gpstr]);
        axis([0 60 0 50]);
        if makeLeg
            legend('expected individual behavior','empirical group behavior','disaster trajectory','95% confidence');
        end
        %suptitle(['Shelter Capacity ' num2str(ss) ', Groups of ' num2str(groupSize)]);
        title(trial_conv(tr,ss,groupSize,groupProtocol));
        hold off; 
end

end

RMSEs = rmse(ts_err');
 % LOOP OVER GROUP SIZE AND protocol?


%% MAKE SUMMARY PLOTS

% for a particular shelter capacity
ss_list = 25;
for ss = ss_list;
figure(find(ss_list==ss));
sizes = [5,25];
sizes = sizes(sizes<=ss);
tbins = 0:1:nts-1;
for groupSize = sizes;
    gix = find(sizes==groupSize);
  for tp=[1,2]
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
sizes = [5,25];
sizes = sizes(sizes<=ss);
load(['ind_grp_crossval_ss' num2str(ss) 't.mat']); % t suffix = crossval; s suffix = sims
load(['ind_grp_crossval_ss' num2str(ss) 's.mat']); % t suffix = crossval; s suffix = sims
for groupSize = sizes;
    gix = find(sizes==groupSize);
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
sizes = [5,25];
sizes = sizes(sizes<=ss);
load(['ind_grp_crossval_ss' num2str(ss) '.mat']);
figure; 
for groupSize = sizes;
    gix = find(sizes==groupSize);
    subplot(1,numel(sizes),gix);
    eval(['plot(tbins,RMSE_ftg' num2str(groupSize) ');']); hold all;
    eval(['plot(tbins,RMSE_ltg' num2str(groupSize) ');']); 
    eval(['plot(tbins,RMSE_mr' num2str(groupSize) ');']); 
    xlabel('time step'); ylabel('RMS prediction error in cum. evacuations');
    legend('FTG','LTG','MV'); axis([0 nts 0 ss]);
    title(['shelter capacity ' num2str(ss) ', group size ' num2str(groupSize)]);
end
