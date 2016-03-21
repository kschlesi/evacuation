function leg = numvec2cell(ib)
% TWO 

if length(ib)~=numel(ib)
    error('must be one-dimensional input');
end

ib = ib(:);
leg = cell(numel(ib),1);
maxL = size(num2str(ib),2);
    for k=1:numel(ib)
        if numel(num2str(ib(k)))<maxL
            leg{k} = [blanks(maxL-numel(num2str(ib(k)))) num2str(ib(k))];
        else
            leg{k} = num2str(ib(k));
        end
    end

end