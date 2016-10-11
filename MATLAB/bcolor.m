function h = bcolor(inmat,varargin)
% provides a balanced color plot (no row/cols left out) with no edge lines
    if ~ismatrix(inmat)
        error('input matrix must be two-dimensional'); 
    end
    
    % find pad value
    pad = mean(mean(inmat));
    if numel(varargin)==2
        varargin{1} = [varargin{1}(:);...
                       varargin{1}(end)+diff(varargin{1}(end-1:end))];
        varargin{2} = [varargin{2}(:);...
                       varargin{2}(end)+diff(varargin{2}(end-1:end))];
    end
    h = pcolor(varargin{:},padarray(inmat,[1 1],pad,'post'));
    set(h, 'EdgeColor', 'none');
    
    % reset axis labels to coincide with centers of appropriate
    % columns/rows
    interp = @(vec_) vec_(1:end-1) + diff(vec_)./2;
    if numel(varargin)==2
        set(gca, 'XTick', interp(varargin{1}), 'XTickLabel', varargin{1}(1:end-1));
        set(gca, 'YTick', interp(varargin{2}), 'YTickLabel', varargin{2}(1:end-1));
    else
        [ysz,xsz] = size(inmat);
        set(gca, 'XTick', interp(1:1:xsz+1), 'XTickLabel', 1:1:xsz);
        set(gca, 'YTick', interp(1:1:ysz+1), 'YTickLabel', 1:1:ysz);
    end
    
end