function [ax,loop,php] = plot_group_ind_samples(cum,ptraj_scaled,Cgrp,tgrp,...
                                                tbins,Cbin,tr,p)
%plot_group_ind_samples(C1(tr,:),Q1(tr,:).*N,Cgrp,tgrp,tbins,Cbin,tr,p,gP)

[ax, loop, php] = plot(cum,'k--'); hold on; % observed in experiment
plot(ptraj_scaled,':');   % Phit trajectory
for i=1:10
	plot(tgrp(:,i),Cgrp(:,i),'-o'); hold all;
	%plot(t_mat(r_mat(:,i)~=0,i),ix(r_mat(r_mat(:,i)~=0,i))-1,':o'); hold all;
end
plot(tbins,mean(Cbin,2),'-k');
title(['trial ' num2str(tr) ', ' num2str(p) ' samples']); axis([0 60 0 50])
legend('empirical data','Phit trajectory (scaled by N)',...
	gP,...%'individual',
	'location','northwest');
xlabel('time'); ylabel('cumulative no. evacuated');

end

