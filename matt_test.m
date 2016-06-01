% matt's idea test...

% cmat = [1 2 2 3 0 1 0;
%            0 3 1 0 1 2 1;
%            0 0 1 0 0 0 1;
%            0 3 2 0 1 2 2;
%            1 2 1 3 0 0 1;
%            0 0 1 0 0 0 1;
%            0 0 0 0 0 0 0];
%        
% cmat = cmat + cmat';       

%% set up test case (in directory \lifespan)
nreg = 11*11+8*7;
pmat = test_case_matt();
cmat = zeros(nreg,nreg);
for sl = 1:size(pmat,4)
    e1 = sub2ind_tcm(pmat(:,1,1,sl),pmat(:,2,1,sl),[1;2]);
    e2 = sub2ind_tcm(pmat(:,1,2,sl),pmat(:,2,2,sl),[1;2]); 
    cmat(sub2ind([nreg,nreg],[e1;e1],[e2;flipud(e2)])) = ...
        cmat(sub2ind([nreg,nreg],[e1;e1],[e2;flipud(e2)])) + 1;
end
cmat = cmat + cmat';
figure; bcolor(cmat);
colorbar;
title('connectivity matrix');

[Ccomm,Q,Chist] = genlouvainREPs(cmat,1,1);
csort = sortentry(cmat,'row',0,Ccomm);
csort = sortentry(csort,'col',0,Ccomm);
figure; bcolor(csort);
colorbar;
title('sorted by communities');

