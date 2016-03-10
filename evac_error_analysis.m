%%

% what are typical differences between neighboring likelihoods?
load('bayesian_pl_ss_ind_fine_err-40-13.mat');
ls = size(A.posteriors); ls(1) = ls(1)-1;
likelihoods = zeros(ls);

% load and plot likelihoods
for i = 1:ls(1)
    likelihoods(i,:,:) = A.posteriors(i+1,:,:)-A.posteriors(i,:,:);
    figure; surf(squeeze(likelihoods(i,:,:))); title(['likelihoods, trial ' num2str(i)]);
end
figure; surf(squeeze(A.posteriors(end,:,:))); title('posterior prob, all trials');

% look at diffs and their distributions
for i=1:ls(1)
    figure(4); hold all;
    diffA=(diff(likelihoods(i,:,:),[],2)); plot(abs(unique(diffA))); hold off;
    disp(min(abs(unique(diffA))));
    figure(5); hold all;
    diffB=(diff(likelihoods(i,:,:),[],3)); plot(abs(unique(diffB))); hold off;
    disp(min(abs(unique(diffB))));
end
