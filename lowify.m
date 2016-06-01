function mat = lowify(mat)

% check that mat is numeric and integers
oldid = unique(mat);
n = numel(oldid);
id = 1:n;
for i=1:n
    mat(mat==oldid(i)) = id(i);
end

end