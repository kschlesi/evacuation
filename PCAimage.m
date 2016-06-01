%% pca compression of image

% original image
matorig = A(1).adj{1};

% data correction
matcorr = bsxfun(@rdivide,bsxfun(@minus,matorig,mean(matorig)),std(matorig));

% covariance matrix
covmat = matcorr'*matcorr;

% eigenvalues and vectors
[evec,lvec] = eig(covmat);
figure(2); hold on;
subplot(1,2,1); bcolor(evec);
subplot(1,2,2); plot(diag(lvec));
suptitle('eig decomp of cov matrix');
hold off;

% determine cutoff dimension & reconstruct
pcts = 0:0.2:1;
cumvars = cumsum(diag(lvec)./sum(diag(lvec)));
for pctvar=pcts
    switch pctvar
        case 0, cutoff = length(cumvars);
        case 1, cutoff = 1;
        otherwise, cutoff = find(cumvars>(1-pctvar),1,'first');
    end
    fvec = evec(:,cutoff:end);
    cdata = matcorr*fvec;
    disp(size(cdata));
    matrec = cdata*fvec';
    matrec = bsxfun(@plus,bsxfun(@times,matrec,std(matorig)),mean(matorig));
    figure(2+find(pcts==pctvar,1,'first')); hold on;
    subplot(1,2,1); bcolor(matorig); 
    title('original matrix');
    subplot(1,2,2); bcolor(matrec);
    title(['reconstructed with pctvar = ' num2str(pctvar) ' (' num2str(size(cdata,2)) ' comp.)']);
    hold off;
end