% ECE594V Final Project
% CSE for clustering in n-1-space

%% first set necessary information (lifespan data -- memory)
inpath = '/Users/kimberly/Documents/lifespan/NNs/';
totalsubjs = 108;
missing_subjs = [46;64;81;82]; % final list
ts_run = 316; nruns = 3;
ts = 52;  % timesteps per window
load('lifespan_subjs.mat');    % get subject IDs ('subj_IDs' 108x1 cell)
t = floor(ts_run/ts)*nruns;
dataset = 'lifespan';

p_plot = {};
% possible member functions of p_plot:
%remove_vis = 0;         % should tagged nodes be removed?
%rm_tag = '_nv';         % sets tag for removing nodes (visual tag = '_nv')

%% or set necessary information (officer data -- rest)
inpath = '/Users/kimberly/Documents/officer_NNs/';
totalsubjs = 86;
missing_subjs = [];
t = 4; ts=0; ts_run=0;
load('officer_subject_names.mat');
subj_IDs = hyperedge_subjs;    % get subject IDs ('subj_IDs' 86x1 cell)
dataset = 'officer';


%% load subject info, adjacency matrices for clustering
ib = loadib(totalsubjs,missing_subjs); % load ib
A = load_adjs(t,ib,subj_IDs,ts,ts_run,inpath,dataset);
n = length(A(1).adj{1});

% hybrid atlas approx locations for BrainView
pos = load('hybrid_centroids.mat');  
pos = squeeze(mean(pos.hybrid_c,2));
pos = squeeze(pos(1,:,:))';


%% convert A.adj to B, the Q-maximization matrix.
gamma = 1;
tmp = cell(size(A));
[A(:).B] = deal(tmp{:});
for k=1:length(A)
  for T=1:t
    ki = sum(A(k).adj{T}); % vector of node degrees
    twom = sum(ki); % m = number of edges in system
    A(k).B{T} = A(k).adj{T} - gamma*(ki'*ki)/twom;
  end
end

%% convert A to E with CSE mapping into n-1 space
usenull = false;
E = A;
tmp = cell(size(E));
[E(:).effdim] = deal(tmp{:});
for k=1:length(E)
    disp(['CSE: subject ' num2str(ib(k))]);
    for T=1:t
        if ~usenull
            [E(k).adj{T}, ed] = CSE(A(k).adj{T});
        else
            [E(k).adj{T}, ed] = CSE(A(k).B{T});
        end
        E(k).effdim = [E(k).effdim; ed];
    end    
end

%% plot effective dimensions
figure;
for k=1:length(E)
    plot(E(k).effdim); hold all;
end
xlabel('run'); ylabel('effective dimension of A');

%% fit gaussian mixture models
kmax = 5; Tmax = 1; Mpick = 5;
Mopt = zeros(kmax,Tmax);
Copt = cell(kmax,Tmax);
for k=1:kmax
  for T=1:Tmax
ed = floor(E(k).effdim(T));
bics = zeros(ed,1);
nQ = zeros(ed,1);
gms = cell(ed,1);
for M=1:ed
    isbic = true;
    try
        gms{M} = fitgmdist(E(k).adj{T}(1:ed,:)',M,'CovarianceType','diagonal');
    catch err
        if strcmp(err.identifier,'stats:gmdistribution:IllCondCovIter')
          disp(['largest M: ' num2str(M-1)]);
          isbic = false;
          if M-1 > 20
              break
          end
        end
    end
    if isbic
        bics(M) = gms{M}.BIC;
        nQ(M) = newmansQ(A(k).adj{T},cluster(gms{M},E(k).adj{T}(1:ed,:)'));
    end
end

bix = find(bics~=0);
bics = bics(bics~=0);
figure(1); hold on;
plot(bix,bics,'.-'); 
hold off;
figure(2); hold on;
plot(bix,nQ(bix),'.-');
hold off;

[~,Mopt(k,T)] = min(bics);
Copt{k,T} = cluster(gms{Mopt(k,T)},E(k).adj{T}(1:ed,:)');
Copt{k,T} = cluster(gms{Mpick},E(k).adj{T}(1:ed,:)');

figure;
csort = sort(Copt{k,T});
Asort = sortentry(sortentry(A(k).adj{T},'row',0,Copt{k,T}),'col',0,Copt{k,T});
firsts = [0;diff(csort)];
Asort(firsts>0,:)=0; % make dividing lines
Asort(:,firsts>0)=0;
bcolor(Asort);

BrainView(pos,Copt{k,T},0);
title(['Optimal CSE GMM clustering, ' num2str(Mopt(k,T)) ' clusters']);

BrainView(pos,cluster(gms{Mpick},E(k).adj{T}(1:ed,:)'),0);
title(['CSE GMM clustering, ' num2str(Mpick) ' clusters']);
  end
end

%% fit k-means clustering models

kmax = 5; Tmax = 1; Mpick = 5;
Mkopt = zeros(kmax,Tmax);
Ckopt = cell(kmax,Tmax);
Cf = cell(kmax,Tmax);
for k=1:kmax
  for T=1:Tmax
edk = floor(E(k).effdim(T));
kms = cell(edk,1);
C = cell(edk,1);
bick = zeros(edk,1);
nQk = zeros(edk,1);
    for M=1:edk
        % fit k-means with M clusters
        [kms{M},C{M}] = kmeans(E(k).adj{T}(1:edk,:)',M);
        % calc bic and save in bick
        %bick(M) = mean(sumd./histc(kms{M},1:M));
        bick(M) = kBIC(E(k).adj{T}(1:edk,:)',kms{M},C{M});
        nQk(M) = newmansQ(A(k).adj{T},kms{M});
    end
    
    bix = find(bick~=0);
    bick = bick(bick~=0);
    figure(1); hold on;
    plot(bix,bick,'.-'); 
    plot(bix,nQk,'.-');hold off;

    [~,Mkopt(k,T)] = min(bick);
    Ckopt{k,T} = kms{Mkopt(k,T)};
    Cf{k,T} = C{Mkopt(k,T)};
    
figure;
csort = sort(Ckopt{k,T});
Asort = sortentry(sortentry(A(k).adj{T},'row',0,Ckopt{k,T}),'col',0,Ckopt{k,T});
firsts = [0;diff(csort)];
Asort(firsts>0,:)=0; % make dividing lines
Asort(:,firsts>0)=0;
bcolor(Asort);

BrainView(pos,Ckopt{k,T},0);
title(['Optimal CSE k-means clustering, ' num2str(Mopt(k,T)) ' clusters']);

BrainView(pos,kms{Mpick},0);
title(['CSE k-means clustering, ' num2str(Mpick) ' clusters']);    
  end
end

%% compute Q for a given C and A



%%



%%
% Load results of Genlouvain clustering instead for comparison

% [gamma, omega] values and timesteps per window
%     goRange = [1.15, 0.005]; ts = 316;  stats_tag = 316;   
%     goRange = [1.15,0.001]; ts = 316; stats_tag = '316alt1';
%     goRange = [1   ,0.005]; ts = 316; stats_tag = '316alt2';
%     goRange = [1   ,0.001]; ts = 316; stats_tag = '316alt3';
goRange = [1.2, 0.05]; ts = 52; stats_tag = '52';