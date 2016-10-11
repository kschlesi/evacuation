function [n,edges,bin] = histcounts(varargin)
%HISTCOUNTS  Histogram Bin Counts.
%   [N,EDGES] = HISTCOUNTS(X) partitions the values in X into bins, and 
%   returns the count in each bin, as well as the bin edges. HISTCOUNTS
%   determines the bin edges using an automatic binning algorithm that 
%   returns uniform bins of a width that is chosen to cover the range of 
%   values in X and reveal the shape of the underlying distribution. 
%
%   N(k) will count the value X(i) if EDGES(k) <= X(i) < EDGES(k+1). The 
%   last bin will also include the right edge such that N(end) will count
%   X(i) if EDGES(end-1) <= X(i) <= EDGES(end).
%
%   [N,EDGES] = HISTCOUNTS(X,M), where M is a scalar, uses M bins.
%
%   [N,EDGES] = HISTCOUNTS(X,EDGES), where EDGES is a vector, specifies the 
%   edges of the bins.
%
%   [N,EDGES] = HISTCOUNTS(...,'BinWidth',BW) uses bins of width BW. To 
%   prevent from accidentally creating too many bins, a limit of 65536 bins
%   can be created when specifying 'BinWidth'. If BW is too small such that 
%   more than 65536 bins are needed, HISTCOUNTS uses wider bins instead.
%
%   [N,EDGES] = HISTCOUNTS(...,'BinLimits',[BMIN,BMAX]) bins only elements 
%   in X between BMIN and BMAX inclusive, X(X>=BMIN & X<=BMAX).
%
%   [N,EDGES] = HISTCOUNTS(...,'Normalization',NM) specifies the
%   normalization scheme of the histogram values returned in N. NM can be:
%                  'count'   Each N value is the number of observations in 
%                            each bin, and SUM(N) is equal to NUMEL(X).
%                            This is the default.
%            'probability'   Each N value is the relative number of 
%                            observations (number of observations in bin / 
%                            total number of observations), and SUM(N) is  
%                            equal to 1.
%           'countdensity'   Each N value is the number of observations in 
%                            each bin divided by the width of the bin. 
%                    'pdf'   Probability density function estimate. Each N 
%                            value is, (number of observations in bin) / 
%                            (total number of observations * width of bin).
%               'cumcount'   Each N value is the cumulative number of 
%                            observations in each bin and all previous bins. 
%                            N(end) is equal to NUMEL(X).
%                    'cdf'   Cumulative density function estimate. Each N 
%                            value is the cumulative relative number of 
%                            observations in each bin and all previous bins. 
%                            N(end) is equal to 1.
%
%   [N,EDGES] = HISTCOUNTS(...,'BinMethod',BM), uses the specified automatic 
%   binning algorithm to determine the number and width of the bins.  BM can be:
%                   'auto'   The default 'auto' algorithm chooses a bin 
%                            width to cover the data range and reveal the 
%                            shape of the underlying distribution.
%                  'scott'   Scott's rule is optimal if X is close to being 
%                            normally distributed, but is also appropriate
%                            for most other distributions. It uses a bin width
%                            of 3.49*STD(X(:))*NUMEL(X)^(-1/3).
%                     'fd'   The Freedman-Diaconis rule is less sensitive to 
%                            outliers in the data, and may be more suitable 
%                            for data with heavy-tailed distributions. It 
%                            uses a bin width of 2*IQR(X(:))*NUMEL(X)^(-1/3), 
%                            where IQR is the interquartile range.
%               'integers'   The integer rule is useful with integer data, 
%                            as it creates a bin for each integer. It uses 
%                            a bin width of 1 and places bin edges halfway 
%                            between integers. To prevent from accidentally 
%                            creating too many bins, a limit of 65536 bins 
%                            can be created with this rule. If the data 
%                            range is greater than 65536, then wider bins 
%                            are used instead.
%                'sturges'   Sturges' rule is a simple rule that is popular
%                            due to its simplicity. It chooses the number of
%                            bins to be CEIL(1 + LOG2(NUMEL(X))).
%                   'sqrt'   The Square Root rule is another simple rule 
%                            widely used in other software packages. It 
%                            chooses the number of bins to be 
%                            CEIL(SQRT(NUMEL(X))).
%
%   [N,EDGES,BIN] = HISTCOUNTS(...) also returns an index array BIN, using 
%   any of the previous syntaxes. BIN is an array of the same size as X 
%   whose elements are the bin indices for the corresponding elements in X. 
%   The number of elements in the kth bin is NNZ(BIN==k), which is the same 
%   as N(k). A value of 0 in BIN indicates an element which does not belong 
%   to any of the bins (for example, a NaN value).
%
%   Class support for inputs X, EDGES:
%      float: double, single
%      integers: uint8, int8, uint16, int16, uint32, int32, uint64, int64
%
%   See also HISTOGRAM

%   Copyright 1984-2014 The MathWorks, Inc.

narginchk(1,inf);
x = varargin{1};

validateattributes(x,{'numeric'},{'real'}, mfilename, 'x', 1)

opts = parseinput(varargin(2:end));

if isempty(opts.BinLimits)  % Bin Limits is not specified
    if ~isempty(opts.BinEdges)
        xc = x(x>=opts.BinEdges(1) & x<=opts.BinEdges(end));
        edges = opts.BinEdges;
    else
        if isinteger(x)
            % for integers, the edges are doubles
            xc = x(:);
            minx = double(min(xc));
            maxx = double(max(xc));
        else
            xc = x(isfinite(x));
            minx = min(xc);  % exclude Inf and NaN
            maxx = max(xc);
        end
        xrange = maxx - minx;
        if ~isempty(opts.NumBins)
            numbins = double(opts.NumBins);
            edges = binpicker(minx,maxx,numbins,xrange/numbins);
        elseif ~isempty(opts.BinWidth)
            if ~isempty(minx)
                % Do not create more than maximum bins.
                MaximumBins = getmaxnumbins();
                binWidth = max(opts.BinWidth, xrange/MaximumBins);
                leftEdge = binWidth*floor(minx/binWidth);
                nbins = max(1,ceil((maxx-leftEdge) ./ binWidth));
                edges = leftEdge + (0:nbins) .* binWidth; % get exact multiples
            else
                edges = cast([0 opts.BinWidth], 'like', xrange);
            end
        else    % BinMethod specified
            if strcmp(opts.BinMethod, 'auto')
                if ~isempty(xc) && (isinteger(xc) || isequal(round(xc),xc))...
                        && xrange <= 50 && maxx <= flintmax(class(maxx))/2 ...
                        && minx >= -flintmax(class(minx))/2
                    opts.BinMethod = 'integers';
                else
                    opts.BinMethod = 'scott';
                end
            end
            switch opts.BinMethod
                case 'scott'
                    binwidth = scottsrule(xc);
                    edges = binpicker(minx,maxx,[],binwidth);
                case 'fd'
                    binwidth = fdrule(xc);
                    edges = binpicker(minx,maxx,[],binwidth);
                case 'integers'
                    if ~isempty(minx)
                        if maxx > flintmax(class(maxx))/2 || ...
                                minx < -flintmax(class(minx))/2
                            error(message('MATLAB:histcounts:InputOutOfIntRange'));
                        end
                        binwidth = integerrule(xc);
                        minx = binwidth*round(minx/binwidth); % make the edges bin width multiples
                        maxx = binwidth*round(maxx/binwidth);
                        edges = (floor(minx)-.5*binwidth):binwidth:(ceil(maxx)+.5*binwidth);
                    else
                        edges = cast([-0.5 0.5], 'like', xrange);
                    end
                case 'sqrt'
                    nbins = sqrtrule(xc);
                    binwidth = (maxx-minx)/nbins;
                    edges = binpicker(minx,maxx,[],binwidth);
                case 'sturges'
                    binwidth = xrange/sturgesrule(xc);
                    edges = binpicker(minx,maxx,[],binwidth);
                    
            end
        end
    end

else   % BinLimits specified
    if isinteger(opts.BinLimits)
        % for integers, the edges are doubles
        minx = double(opts.BinLimits(1));
        maxx = double(opts.BinLimits(2));
    else
        minx = opts.BinLimits(1);  
        maxx = opts.BinLimits(2);
    end
    xrange = maxx - minx;
    xc = x(x>=minx & x<=maxx);
    if ~isempty(opts.NumBins)
        numbins = double(opts.NumBins); 
        edges = [minx + (0:numbins-1).*(xrange/numbins), maxx];            
    elseif ~isempty(opts.BinWidth)
        % Do not create more than maximum bins.
        MaximumBins = getmaxnumbins();
        binWidth = max(opts.BinWidth, xrange/MaximumBins);
        edges = minx:binWidth:maxx; 
        if edges(end) < maxx || isscalar(edges)
            edges = [edges maxx];
        end

    else    % BinMethod specified
        if ~isempty(xc)
            if strcmp(opts.BinMethod, 'auto')
                if ~isempty(xc) && (isinteger(xc) || isequal(round(xc),xc)) ...
                        && xrange <= 50 && maxx <= flintmax(class(maxx))/2 ...
                        && minx >= -flintmax(class(minx))/2
                    opts.BinMethod = 'integers';
                else
                    opts.BinMethod = 'scott';
                end
            end
            switch opts.BinMethod
                case 'scott'
                    binwidth = scottsrule(xc);
                    nbins = ceil((maxx-minx)/binwidth);
                    edges = linspace(minx,maxx,nbins);
                case 'fd'
                    binwidth = fdrule(xc);
                    nbins = ceil((maxx-minx)/binwidth);
                    edges = linspace(minx,maxx,nbins);
                case 'integers'
                    if maxx > flintmax(class(maxx))/2 || ...
                            minx < -flintmax(class(minx))/2
                        error(message('MATLAB:histcounts:InputOutOfIntRange'));
                    end
                    binwidth = integerrule(xc);
                    minxi = binwidth*ceil(minx/binwidth)+0.5;
                    maxxi = binwidth*floor(maxx/binwidth)-0.5;
                    edges = [minx minxi:binwidth:maxxi maxx]; 
                case 'sqrt'
                    nbins = sqrtrule(xc);
                    edges = linspace(minx,maxx,nbins);
                case 'sturges'
                    nbins = sturgesrule(xc);
                    edges = linspace(minx,maxx,nbins);
                    
            end
        else
            if strcmp(opts.BinMethod, 'integers')
                minxi = ceil(minx)+0.5;
                maxxi = floor(maxx)-0.5;
                edges = [minx minxi:maxxi maxx];
            else
                edges = [minx maxx];
            end
        end
       
    end
end

edges = full(edges); % make sure edges are non-sparse
if nargout <= 2
    n = histcountsmex(x,edges);
else
    [n,bin] = histcountsmex(x,edges);
end

switch opts.Normalization
    case 'countdensity'
        n = n./double(diff(edges));
    case 'cumcount'
        n = cumsum(n);
    case 'probability'
        n = n / numel(xc);
    case 'pdf'
        n = n/numel(xc)./double(diff(edges));
    case 'cdf'
        n = cumsum(n / numel(xc));
end
    
end

function opts = parseinput(input)

opts = struct('NumBins',[],'BinEdges',[],'BinLimits',[],...
    'BinWidth',[],'Normalization','count','BinMethod','auto');
funcname = mfilename;

% Parse second input in the function call
if ~isempty(input)
    in = input{1};
    numposinput = 0;
    if isnumeric(in)
        if isscalar(in)
            validateattributes(in,{'numeric'},{'integer', 'positive'}, ...
                funcname, 'm', 2)
            opts.NumBins = in;   
            opts.BinMethod = [];
        else
            validateattributes(in,{'numeric'},{'vector','nonempty', ...
                'real', 'nondecreasing'}, funcname, 'edges', 2)
            opts.BinEdges = in;
            opts.BinMethod = [];
        end
        input(1) = [];
        numposinput = 1;
    end
    
    % All the rest are name-value pairs
    inputlen = length(input);
    if rem(inputlen,2) ~= 0
        error(message('MATLAB:histcounts:ArgNameValueMismatch'))
    end
    
    for i = 1:2:inputlen
        name = validatestring(input{i}, {'NumBins', 'BinEdges', 'BinWidth', 'BinLimits', ...
            'Normalization', 'BinMethod'}, i+1+numposinput);
        
        value = input{i+1};
        switch name
            case 'NumBins'
                validateattributes(value,{'numeric'},{'scalar', 'integer', ...
                    'positive'}, funcname, 'NumBins', i+2+numposinput)
                opts.NumBins = value;
                if ~isempty(opts.BinEdges)
                    error(message('MATLAB:histcounts:InvalidMixedBinInputs'))
                end
                opts.BinMethod = [];
                opts.BinWidth = [];
            case 'BinEdges'
                validateattributes(value,{'numeric'},{'vector','nonempty', ...
                    'real', 'nondecreasing'}, funcname, 'BinEdges', i+2+numposinput);
                opts.BinEdges = reshape(value,1,[]);
                opts.BinMethod = [];
                opts.NumBins = [];
                opts.BinWidth = [];
                opts.BinLimits = [];
            case 'BinWidth'
                validateattributes(value, {'numeric'}, {'scalar', 'real', ...
                    'positive', 'finite'}, funcname, 'BinWidth', i+2+numposinput);
                opts.BinWidth = value;
                if ~isempty(opts.BinEdges)
                    error(message('MATLAB:histcounts:InvalidMixedBinInputs'))
                end
                opts.BinMethod = [];
                opts.NumBins = [];
            case 'BinLimits'
                validateattributes(value, {'numeric'}, {'numel', 2, 'real', ...
                    'nondecreasing', 'finite'}, funcname, 'BinLimits', i+2+numposinput)
                opts.BinLimits = value;
                if ~isempty(opts.BinEdges)
                    error(message('MATLAB:histcounts:InvalidMixedBinInputs'))
                end
            case 'Normalization'
                opts.Normalization = validatestring(value, {'count', 'countdensity', 'cumcount',...
                    'probability', 'pdf', 'cdf'}, funcname, 'Normalization', i+2+numposinput);
            otherwise % 'BinMethod'
                opts.BinMethod = validatestring(value, {'auto','scott', 'fd', ...
                    'integers', 'sturges', 'sqrt'}, funcname, 'BinMethod', i+2+numposinput);
                if ~isempty(opts.BinEdges)
                    error(message('MATLAB:histcounts:InvalidMixedBinInputs'))
                end
                opts.BinWidth = [];
                opts.NumBins = [];
        end
    end
end  
end

function mb = getmaxnumbins
mb = 65536;  %2^16
end

function binwidth = scottsrule(x)
% Scott's normal reference rule
if isinteger(x)
    x = double(x);
else
    x = x(isfinite(x));
end
binwidth = 3.5*std(x)/(numel(x)^(1/3));
end

function y = iqr(x)
% IQR Compute interquartile range.
n = numel(x);
F = ((1:n)'-.5) / n;
if isinteger(x)
    x = double(x);
end
y = diff(interp1q(F, sort(x(:)), [.25; .75]));
end

function binwidth = fdrule(x)
n = numel(x);
xrange = max(x(:)) - min(x(:));
if n > 1
    % Guard against too small an IQR.  This may be because there
    % are some extreme outliers.
    iq = max(iqr(x),double(xrange)/10);
    binwidth = 2 * iq * n^(-1/3);
else
    binwidth = 1;
end
end

function binwidth = integerrule(x)
MaximumBins = getmaxnumbins();
xscale = max(abs(x(:)));
xrange = max(x(:)) - min(x(:));
if xrange > MaximumBins
    % If there'd be more than maximum bins, center them on an appropriate
    % power of 10 instead.
    binwidth = 10^ceil(log10(xrange/MaximumBins));
elseif isfloat(x) && eps(xscale) > 1
    % If a bin width of 1 is effectively zero relative to the magnitude of
    % the endpoints, use a bigger power of 10.
    binwidth = 10^ceil(log10(eps(xscale)));
else
    % Otherwise bins are centered on integers.
    binwidth = 1;
end
end

function nbins = sturgesrule(x)
nbins = max(ceil(log2(numel(x))+1),1);
end

function nbins = sqrtrule(x)
nbins = max(ceil(sqrt(numel(x))),1);
end

function edges = binpicker(xmin,xmax,nbins,rawBinWidth)
% BINPICKER Choose histogram bins.

if ~isempty(xmin)
    xscale = max(abs([xmin,xmax]));
    xrange = xmax - xmin;
    
    % Make sure the bin width is not effectively zero.  Otherwise it will never
    % amount to anything, which is what we knew all along.
    rawBinWidth = max(rawBinWidth, eps(xscale));
    
    % If the data are not constant, place the bins at "nice" locations
    if xrange > max(sqrt(eps(xscale)), realmin(class(xscale)))
        % Choose the bin width as a "nice" value.
        powOfTen = 10.^floor(log10(rawBinWidth)); % next lower power of 10
        relSize = rawBinWidth / powOfTen; % guaranteed in [1, 10)
        
        % Automatic rule specified
        if isempty(nbins)
        if  relSize < 1.5
            binWidth = 1*powOfTen;
        elseif relSize < 2.5
            binWidth = 2*powOfTen;
        elseif relSize < 4
            binWidth = 3*powOfTen;
        elseif relSize < 7.5
            binWidth = 5*powOfTen;
        else
            binWidth = 10*powOfTen;
        end
        
            % Put the bin edges at multiples of the bin width, covering x.  The
            % actual number of bins used may not be exactly equal to the requested
            % rule. 
            leftEdge = min(binWidth*floor(xmin ./ binWidth), xmin);
            nbinsActual = max(1, ceil((xmax-leftEdge) ./ binWidth));
            rightEdge = max(leftEdge + nbinsActual.*binWidth, xmax);
            
        % Number of bins specified
        else            
            % temporarily set a raw binWidth to a nice power of 10. 
            % binWidth will be set again to a different value if more than 
            % 1 bin.
            binWidth = powOfTen*floor(relSize);
            % Set the left edge at multiples of the raw bin width.
            % Then adjust bin width such that all bins are of the same
            % size and xmax fall into the rightmost bin.
            leftEdge = min(binWidth*floor(xmin ./ binWidth), xmin);
            if nbins > 1
                ll = (xmax-leftEdge)/nbins;  % binWidth lower bound, xmax
                                             % on right edge of last bin
                ul = (xmax-leftEdge)/(nbins-1);  % binWidth upper bound,
                                          % xmax on left edge of last bin
                p10 = 10^floor(log10(ul-ll));
                binWidth = p10*ceil(ll./p10);  % binWidth-ll < p10 <= ul-ll
                                               % Thus, binWidth < ul                 
            end
            
            nbinsActual = nbins;
            rightEdge = max(leftEdge + nbinsActual.*binWidth, xmax);
        end
        
    else % the data are nearly constant
        % For automatic rules, use a single bin.
        if isempty(nbins)
            nbins = 1;
        end
        
        % There's no way to know what scale the caller has in mind, just create
        % something simple that covers the data.

        % Make the bins cover a unit width, or as small an integer width as
        % possible without the individual bin width being zero relative to
        % xscale.  Put the left edge on an integer or half integer below
        % xmin, with the data in the middle 50% of the bin.  Put the right
        % edge similarly above xmax.
        binRange = max(1, ceil(nbins*eps(xscale)));
        leftEdge = floor(2*(xmin-binRange./4))/2;
        rightEdge = ceil(2*(xmax+binRange./4))/2;

        binWidth = (rightEdge - leftEdge) ./ nbins;
        nbinsActual = nbins;
    end
    
    edges = [leftEdge + (0:nbinsActual-1).*binWidth, rightEdge];
else
    % empty input
    if ~isempty(nbins)
        edges = cast(0:nbins, 'like', xmin);
    else
        edges = cast([0 1], 'like', xmin);
    end
end

end