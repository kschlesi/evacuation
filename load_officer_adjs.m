function A = load_officer_adjs(t,ib,subj_IDs,ts,run_ts,inpath,dataset)

if nargin<7
    dataset = 'lifespan';
end

A = struct('number',cell(numel(ib),1),'ID',cell(numel(ib),1),...
           'adj',cell(numel(ib),1));

switch dataset
    case 'lifespan',
        TperR = floor(run_ts/ts);
        for k=1:numel(ib)
          disp(['subject ' num2str(ib(k))]);
          tempcell = cell(t,1);
          for T = 1:t
            r = floor((T-1)/TperR)+1;
            i = mod(T-1,TperR)+1;
            temp1 = csvread([inpath 'ageNN' num2str(ib(k)) '_run' num2str(r) ...
                         '-' num2str(i) '_tw' num2str(ts) '.txt']);
            temp1(isnan(temp1)) = 0;
            tempcell{T} = temp1;
          end
          A(k).number = ib(k);
          A(k).ID = subj_IDs{ib(k)};
          A(k).adj = tempcell;
        end
    case 'officer',
        for k=1:numel(ib)
          disp(['subject ' num2str(ib(k))]);
          tempcell = cell(t,1);
          for T = 1:t
            
            temp1 = csvread([inpath 'ageNN' num2str(ib(k)) '_run' num2str(r) ...
                         '-' num2str(i) '_tw' num2str(ts) '.txt']);
            temp1(isnan(temp1)) = 0;
            tempcell{T} = temp1;
          end
          A(k).number = ib(k);
          A(k).ID = subj_IDs{ib(k)};
          A(k).adj = tempcell;
        end
        
end

end