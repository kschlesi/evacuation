function makeLeg = fig_format_sum(tr,horm,tr_hit,tr_miss,gPs,groupProtocol)
% creates proper subplot call for 'sum' (summary) figure format
% returns makeLeg (boolean indicator for legend in given subplot)

         assert(find([tr_hit,tr_miss]==tr)==(2-horm));
         makeLeg = false;
         subplot(3,2,2*find(ismemvar(gPs,groupProtocol))-horm);
         if horm && strcmp(groupProtocol,'mR')
            makeLeg = true;
         end
         
end