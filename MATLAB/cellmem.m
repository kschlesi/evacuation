function [comp,xinbs] = cellmem(A,B)
% this takes two cell arrays, A and B. returns a logical vector C the same
% size as A, with C(i) = 1 if the contents of A(i) are identical to the
% contents of ANY cell in B.

comp = zeros(size(A));
xinbs = zeros(size(A));
for a=1:numel(A)
    xinb = zeros(size(B));
    for b=1:numel(B)
        % if either cell (A or B) CONTAINS a char, use strcmp;
        % else use all(equals)
        if ischar(A{a}) || ischar(B{b})
            xinb(b) = strcmp(A{a},B{b});
        else
            try AeqB = A{a}==B{b};
            catch
                AeqB = 0;
            end
            xinb(b) = (all(AeqB(:)) && numel(A{a})==numel(B{b}));
        end
    end
    comp(a) = any(xinb(:));
    if comp(a);
        xinbs(a) = find(~~xinb,1,'first');
    end
end
        
end