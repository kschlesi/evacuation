function [proj, u] = vnique(C,dim)
% takes a multidimensional array C and returns the number of
% unique elements projected along any dimension.

% currently supports only 2D arrays.
u = unique(C);
proj = zeros(size(C,dim),1);
for i=u'
    iproj = ~~sum((C==i),-1*(dim-3));
    proj = proj + iproj(:);
end