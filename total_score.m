function [tot,all1] = total_score(q_model,pvec,Phit,trials,Q1,lossmat)
% computes total score over all trials in 'trials'
all1 = zeros(numel(trials),1);
for i=trials(:)'
    all1(trials==i) = exp_score(q_model,pvec,Phit,Q1(i,:),lossmat);
    if isnan(all1(trials==i))
        disp(['trial ' num2str(i)]);
    end
end
tot = sum(all1);

end