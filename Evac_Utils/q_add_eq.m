function [q_neq,q_eq] = q_add_eq(pvec,q_con)
% need nargin(qcon) = 1 (pvec)
q_neq = q_con(pvec);
q_eq = 0;