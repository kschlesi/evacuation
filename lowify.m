function mat = lowify(mat,ignoreval)

% check that mat is numeric and integers
oldid = unique(mat);
if nargin>1
    oldid = removeval(oldid,ignoreval);
end
n = numel(oldid);
id = 1:n;
for i=1:n
    mat(mat==oldid(i)) = id(i);
end

end