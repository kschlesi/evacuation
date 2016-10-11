function makeLeg = fig_format_many(tr,test_trials,ss,groupProtocol,groupSize)
% creates proper subplot call for 'many' figure format
% returns makeLeg (boolean indicator for legend in given subplot)

          if groupSize==5 && ss==50
            makeLeg = false;
            switch groupProtocol
                case 'fTG',
                  subplot(3,6,find(test_trials==tr));
                case 'mR'
                  subplot(3,6,[11,12,17,18]);
                  makeLeg = true;
                case 'lTG'
                  if find(test_trials==tr)<2  
                    subplot(3,6,find(test_trials==tr)+9);
                  else
                    subplot(3,6,find(test_trials==tr)+11);  
                  end
            end
          end
          if groupSize==25 && ss==50
             makeLeg = false;
             switch groupProtocol
                case 'fTG',
                  subplot(4,5,find(test_trials==tr));
                case 'mR'
                  if find(test_trials==tr)==1  
                  subplot(4,5,5);
                  else
                      if find(test_trials==tr)<7
                          subplot(4,5,find(test_trials==tr)+4);
                      else
                          subplot(4,5,[14,15,19,20]);
                          makeLeg = true;
                      end
                  end
                case 'lTG'
                  if find(test_trials==tr)<4
                    subplot(4,5,find(test_trials==tr)+10);
                  else
                    subplot(4,5,find(test_trials==tr)+12);
                  end
            end  
          end
          if groupSize==5 && ss==25
            makeLeg = false;
            switch groupProtocol
                case 'fTG',
                  subplot(3,7,find(test_trials==tr));
                case 'mR'
                  if find(test_trials==tr)<6
                    subplot(3,7,find(test_trials==tr)+7);
                  else
                    makeLeg = true;
                    subplot(3,7,[13,14,20,21]);
                  end
                case 'lTG'
                  subplot(3,7,find(test_trials==tr)+14);
            end
          end
          if groupSize==25 && ss==25
             makeLeg = false;
             switch groupProtocol
                case 'fTG',
                  subplot(4,5,find(test_trials==tr));
                case 'mR'
                  if find(test_trials==tr)<3  
                    subplot(4,5,find(test_trials==tr)+3);
                  else
                    subplot(4,5,find(test_trials==tr)+3);  
                  end
                case 'lTG'
                  if find(test_trials==tr)<4  
                    subplot(4,5,find(test_trials==tr)+10);
                  else
                    if find(test_trials==tr)<7  
                      subplot(4,5,find(test_trials==tr)+12);
                    else
                      makeLeg = true;
                      subplot(4,5,[14,15,19,20]);  
                    end
                  end
            end  
          end
          
end