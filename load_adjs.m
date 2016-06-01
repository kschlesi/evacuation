function A = load_adjs(t,ib,subj_IDs,ts,run_ts,inpath,dataset)

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
          if t==4
            tasks = [{'rest'},{'attn'},{'word'},{'face'}];
            tasknos = repmat({[]},t,1);
            path1 = [inpath 'N'];
            path2 = ['_' num2str(ib(k)) '.dat'];
          end
          if t==41
            tasks = [repmat({'rest'},3,1);repmat({'att'},12,1);repmat({'word'},13,1);repmat({'face'},13,1)];
            tasknos = num2cell([1:3,1:12,1:13,1:13]');
            path1 = [inpath 'ocd' num2str(ib(k)) '_'];
            path2 = '_tw40_atl194.txt';
          end
          for T = 1:t
            temp1 = csvread([path1,char(tasks(T)),num2str(tasknos{T}),path2]);
            temp1(isnan(temp1)) = 0;
            tempcell{T} = temp1;
          end
          A(k).number = ib(k);
          A(k).ID = subj_IDs{k};
          A(k).adj = tempcell;
        end
        
end

end