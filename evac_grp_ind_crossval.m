%% CROSS-VALIDATION 

%% set training and holdout trials; set shelter space

ss = 50;
train_trials = trial_ix('ind',ss,1,missing);

% for LOOCV on ind trials:
%loo_ix = 3;
%test_trials = train_trials(loo_ix);
%train_trials = removeval(train_trials,test_trials);

% for cross-validation on group trials:
groupProtocol = 'lTG';
groupSize = 25;
test_trials = trial_ix(groupProtocol,ss,groupSize,missing);

makeFigs = true;

%% load necessary data

bins = -0.05:0.1:1.05;  % this gives the bin EDGES for Phit values
[Htrain,Jtrain,Theta,AvgCumEvac,missing,Q1,T1,P1,C1] = load_evac_data(0,train_trials,bins);
[Htest,Jtest] = load_evac_data(0,test_trials,bins);
[ntrials,nts] = size(Q1);
[N,~] = size(P1);

%% first, perform (individual) fit on training trials
% with maximum likelihood estimation

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
    Cix(tbins==ix) = Cix(find(tbins==ix)-1);
  end
end
tot_err(test_trials==tr) = Cexp(end) - C1(tr,end);
ts_err(test_trials==tr,:) = Cexp(Cix)' - C1(tr,:);

if makeFigs
figure; bcolor(P',T,0:1:N); colorbar; hold all; % probability distribution
        plot(T,Cexp);                           % expected value
        plot(tbins,C1(tr,:));                         % empirical value
        plot(tbins,Q1(tr,:).*N,'--');                 % Phit value
        xlabel('time'); ylabel('# evacuated'); title(['trial ' num2str(tr)]);
        hold off;
end

end

%% comparison:

RMSE = sqrt(mean(ts_err.^2,1));

%% MAKE PLOTS

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
