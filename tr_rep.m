function [tr_hit,tr_miss] = tr_rep(ss,gs,gp)

if ss==50
   if gs==5 
       switch gp
           case 'fTG', tr_hit = 122; 
                       tr_miss = 35;
           case 'lTG', tr_hit = 53; 
                       tr_miss = 37;
           case 'mR', tr_hit = 56; 
                      tr_miss = NaN;
       end
   end
   if gs==25
       switch gp
           case 'fTG', tr_hit = 117; 
                       tr_miss = 38;
           case 'lTG', tr_hit = 111; 
                       tr_miss = 77;
           case 'mR', tr_hit = 81; 
                      tr_miss = 106;
       end
   end
end

if ss==25
   if gs==5 
       switch gp
           case 'fTG', tr_hit = 46; 
                       tr_miss = 69;
           case 'lTG', tr_hit = 58; 
                       tr_miss = 119;
           case 'mR', tr_hit = 83; 
                      tr_miss = 85;
       end
   end
   if gs==25
       switch gp
           case 'fTG', tr_hit = 92; 
                       tr_miss = 142;
           case 'lTG', tr_hit = 140; 
                       tr_miss = 80;
           case 'mR', tr_hit = 21; 
                      tr_miss = 33;
       end
   end
end


end