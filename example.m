% This script runs the GBHC algorithm:
% http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0075748
% 
% Bayesian Hierarchical Clustering for Studying Cancer Gene Expression Data with Unknown Statistics
% Korsuk Sirinukunwattana, Richard S. Savage, Muhammad F. Bari, David R. J. Snead,
% Nasir M. Rajpoot

% Update Log :

% 2014_1_14 :
% Include an option of using parallel computing toolbox. 
% Thanks Jonathan Erick <ericksonj@wlu.edu> for pointing this out.

% 2014_1_17 : 
% the second part of the code 'Make some pretty plot' is
% kindly provided by Jonathan Erickson <ericksonj@wlu.edu>. It uses 
% plot_gaussian_ellipsoid.m, from  
% http://www.mathworks.com/matlabcentral/fileexchange/16543-plotgaussianellipsoid
% written by Gautam Vallabha, to make a guassian plot 


%% Perform Clustering 

% Generate Data
R1 = normrnd(0,0.1,[20,2]);
R2 = normrnd(0.5,0.1,[20,2]);
R3 = normrnd(1,0.2,[20,2]);
data = [R1;R2;R3];

% Speed up the algorithms by parallel computing
% Open matlab workers
v = ver;
toolboxName = 'Parallel Computing Toolbox';
PCTavailable = any(strcmp(toolboxName, {v.Name}));
if PCTavailable
    s = matlabpool('size');
    if s == 0
        matlabpool open
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GBHC with hyperparameter optimization globally over the whole tree (TREE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = 10; % number of random starting points for hyperparameters optimisation
[idx1,dendro] = GBHC_TREE(data,m);
% plot dendrogram
probeID = {'1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1',...
           '2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2',...
           '3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3'};
orient = 'left';
figure(1), BHCDendrogram(dendro,0,'labels',probeID,'orientation',orient);

title('Dendrogram Produced by GBHC-TREE')
% write cluster's index
fprintf('\n\n\n Cluster membership inferred by GBHC-TREE:\n')
idx1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GBHC with hyperparameter optimization at every merger (NODE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[idx2,dendro] = GBHC_NODE(data);
% plot dendrogram
probeID = {'1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1',...
           '2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2',...
           '3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3'};
orient = 'left';
figure(2), BHCDendrogram(dendro,0,'labels',probeID,'orientation',orient);
title('Dendrogram Produced by GBHC-NODE')
% write cluster's index
fprintf('\n\n\n Cluster membership inferred by GBHC-NODE:\n')
idx2

% close matlab workers
if PCTavailable
  matlabpool close
end


%% Make some pretty plots.
% 1. plot the discovered data clusters
% GBHC-TREE
figure; subplot(2,2,1)
hold on;
title('Discovered data clusters GBHC-TREE')

jj = find(idx1==1);
plot(data(jj,1), data(jj,2), 'b*')

jj = find(idx1==2);
plot(data(jj,1), data(jj,2), 'r*')

jj = find(idx1==3);
plot(data(jj,1), data(jj,2), 'k*')

% GBHC-NODE
subplot(2,2,2)
hold on;
title('Discovered data clusters GBHC-NODE')

jj = find(idx2==1);
plot(data(jj,1), data(jj,2), 'b*')

jj = find(idx2==2);
plot(data(jj,1), data(jj,2), 'r*')

jj = find(idx2==3);
plot(data(jj,1), data(jj,2), 'k*')

% 2. plot the actual clustered

subplot(2,2,3)
hold on;
title('Actual data clusters')
plot(data(1:20,1), data(1:20,2), 'b*')
plot(data(21:40,1), data(21:40,2), 'r*')
plot(data(41:60,1), data(41:60,2), 'k*')

% Draw Gaussian ellipsoids to show original
% generating components.

h3 = plot_gaussian_ellipsoid([1 1], [0.2 0; 0 0.2].^2, 1); set(h3, 'color', 'k')
h3 = plot_gaussian_ellipsoid([1 1], [0.2 0; 0 0.2].^2, 2); set(h3, 'color', 'k')
h3 = plot_gaussian_ellipsoid([1 1], [0.2 0; 0 0.2].^2, 0.5); set(h3, 'color', 'k')

h2 = plot_gaussian_ellipsoid([0.5 0.5], [0.1 0; 0 0.1].^2, 1); set(h2, 'color', 'r')
h2 = plot_gaussian_ellipsoid([0.5 0.5], [0.1 0; 0 0.1].^2, 0.5); set(h2, 'color', 'r')
h2 = plot_gaussian_ellipsoid([0.5 0.5], [0.1 0; 0 0.1].^2, 2); set(h2, 'color', 'r')

h1 = plot_gaussian_ellipsoid([0 0], [0.1 0; 0 0.1].^2, 0.1); set(h1, 'color', 'b')
h1 = plot_gaussian_ellipsoid([0 0], [0.1 0; 0 0.1].^2, 0.5); set(h1, 'color', 'b')
h1 = plot_gaussian_ellipsoid([0 0], [0.1 0; 0 0.1].^2, 2); set(h1, 'color', 'b')
