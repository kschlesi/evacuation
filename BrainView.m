function [] = BrainView(pos,nodedata,edgedata,varargin)
% This function produces a 3D plot of brain data, given the 3D positions of
% brain nodes and data values associated with each node.
% REQUIRED INPUTS: "pos" : an nx3 matrix containing (x,y,z) coordinates for
%                              each of n nodes
%                  "nodedata" : an nx1 vector containing numerical values
%                                   to be represented at each node
%                  "edgedata" : a matrix containing edge data; set to 0
%                                   (since edge plotting is not supported)
%
% The default node property to be determined by 'nodedata' is node color.
% With no other arguments, each node will be colored according to its value
% in "nodedata" (and colors scaled using the default colormap).
% To change the default node property to be determined by "nodedata"
% include the following arguments: 'DefaultNodeProperty' followed by the
% desired property. 
% Example: for "nodedata" to be interpreted as a vector of node sizes, use:
%     BrainView(pos,nodedata,0,'DefaultNodeProperty','size')
% All possible default node properties: 'color' : node color
%                                       'size' : node size
%                                       'border' : node border color
%                                       'stroke' : node border line stroke
%
% To control a node property other than the default node property, include
% the following arguments: a keyword of the form 'Node<Property>' 
% specifying which property to control, followed by a scalar or an
% n-element vector containing the desired values.
% Example: to control the node border strokes with the nx1 vector "bStroke"
%     use: BrainView(pos,nodedata,0,'NodeStroke',bStroke)
% Example: to additionally make all of the node border colors blue, use:
%    BrainView(pos,nodedata,0,'NodeStroke',bStroke,'NodeBorder','b') or
%    BrainView(pos,nodedata,0,'NodeStroke',bStroke,'NodeBorder',[0 0 1])
%
% VIEW: There are two ways to adjust the view with the 'View' keyword.
% (1) Use a built-in setting -- either 'axial' 'coronal' or 'sagittal':
%       BrainView(pos,nodedata,0,'View','axial')
% (2) Set a custom rotation with a 2-element vector containing the xy-plane
%       rotation angle followed by the elevation angle, for example:
%       BrainView(pos,nodedata,0,'View',[90,45])
%       will rotate 90 degrees in the xy-plane and elevate 45 degrees.
% 
% NOTE: this function automatically creates a new figure with the "figure;" 
% command before plotting. To eliminate this (for example, when creating
% sub-figures), include the keyword 'ExtFigureCmd' as one of the arguments.

[n,dim] = size(pos);
if dim~=3; error('must provide 3-dimensional coordinates'); end;
if isscalar(nodedata); nodedata = repmat(nodedata,n,1); end;

% set node property to be determined by nodedata
if any(ismemvar(varargin,'DefaultNodeProperty'))
    defNode = char(varargin(find(ismemvar(varargin,'DefaultNodeProperty'),1,'first')+1));
else
    defNode = 'color';
end
eval(['node' defNode 's = nodedata;']);
% default node color = black
if ~strcmp(defNode,'color'); nodecolors = repmat('k',n,1); end;

% set edge property to be determined by edgedata
if any(ismemvar(varargin,'DefaultEdgeProperty'))
    defEdge = char(varargin(find(ismemvar(varargin,'DefaultEdgeProperty'),1,'first')+1));
else
    defEdge = 'stroke';
end

% set other possible node properties (size, border color, boreder stroke)
% format: if not default; if passed in, read in & convert scalar to vector
%                         if not passed in, set default
if ~strcmp(defNode,'size')
    if any(ismemvar(varargin,'NodeSize'))
        nodesizes = varargin{find(ismemvar(varargin,'NodeSize'),1,'first')+1};
        if isscalar(nodesizes); nodesizes = repmat(nodesizes,n,1); end;
    else
        nodesizes = 100*ones(n,1); % default node size of 100
    end
end
if ~strcmp(defNode,'border')
    if any(ismemvar(varargin,'NodeBorder'))
        nodeborders = varargin{find(ismemvar(varargin,'NodeBorder'),1,'first')+1};
        if isscalar(nodeborders) || all(size(nodeborders)==[1,3])
            nodeborders = repmat(nodeborders,n,1);
        end
    else
        nodeborders =  repmat('k',n,1); % default black node borded
    end
end
if ~strcmp(defNode,'stroke')
    if any(ismemvar(varargin,'NodeStroke'))
        nodestrokes = varargin{find(ismemvar(varargin,'NodeStroke'),1,'first')+1};
        if isscalar(nodestrokes); nodestrokes = repmat(nodestrokes,n,1); end;
    else
        nodestrokes = 1.25.*ones(n,1); % default node border stroke of 1.25
    end
end

% set other possible edge properties (not yet supported)
if edgedata
   warning('edge plotting not yet supported');
   edgedata = 0;
end

if ~any(ismemvar(varargin,'ExtFigureCmd'));
figure;    
end
% create node plot
for i = 1:n;
    scatter3(pos(i,1),pos(i,2),pos(i,3),nodesizes(i),nodecolors(i,:),'fill',...
        'MarkerEdgeColor',nodeborders(i,:),'LineWidth', nodestrokes(i));
    hold all
end

% create edge plot
if edgedata
    % create edge plot; not yet supported
    warning('Plotting of edge data not supported');    
end

% adjust view if desired
if any(ismemvar(varargin,'View'))
    viewin = varargin{find(ismemvar(varargin,'View'),1,'first')+1};
    if ischar(viewin)
        switch viewin
            case 'sagittal', view(90,0); title('Sagittal View');
            case 'coronal', view(0,0); title('Coronal View');
            case 'axial', view(0,90); title('Axial View');
        end
    else
        if numel(viewin)==2; view(viewin(1),viewin(2)); end;
    end
end    

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Lia,Locb] = ismemvar(A,B,varargin)
% this function works like 'ismember,' but allows operation on cell arrays 
% with both string and non-string cells. 
% *****
% NOTE: if a combination of cell array + non-cell array is passed in that 
% would cause an error in 'ismember' (i.e., the cell array is not a cell 
% array of strings OR the non-cell array is not a string), the DEFAULT
% behavior of this function is to check whether the contents of each cell 
% in the cell array are identical to the ENTIRE non-cell array input.
% EXAMPLES: ismemvar({'one'; 'two'; 1; 2}, 'one') --> [1; 0; 0; 0]
%           ismemvar({'one'; 'two'; [1 2]}, [1 2]) --> [0; 0; 1]
%           ismemvar({'one'; 'two'; 1; 2}, [1 2]) --> [0; 0; 0; 0]

% cases in which it directly mimics 'ismember'
if (~iscell(A) && ~iscell(B)) || ...
    ((ischar(A) || iscellstr(A)) && (ischar(B) || iscellstr(B)))
    [Lia,Locb] = ismember(A,B,varargin{:});
else
    % otherwise.... use cellmem. Locb is NOT supported here.
    if ~iscell(A); A = {A}; end;
    if ~iscell(B); B = {B}; end;
    Lia = cellmem(A,B);
    Locb = [];
    if numel(varargin)
      warning('rows or legacy behavior not supported for these arguments'); 
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function comp = cellmem(A,B)
% this takes two cell arrays, A and B. returns a logical vector C the same
% size as A, with C(i) = 1 if the contents of A(i) are identical to the
% contents of ANY cell in B.

comp = zeros(size(A));
for a=1:numel(A)
    xinb = zeros(size(B));
    for b=1:numel(B)
        % if either cell (A or B) CONTAINS a char, use strcmp;
        % else use all(equals)
        if ischar(A{a}) || ischar(B{b})
            xinb(b) = strcmp(A{a},B{b});
        else
            try AeqB = A{a}==B{b};
            catch oops
                if oops; AeqB = 0; end;
            end
            xinb(b) = (all(AeqB(:)) && numel(A{a})==numel(B{b}));
        end
    end
    comp(a) = any(xinb(:));
end
        
end