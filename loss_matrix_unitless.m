function lossmat = loss_matrix_unitless(rho_a,rho_b,span,intercept)

% default intercept = 0
if nargin<4
    intercept = 0;
end

% default span = 10
if nargin<3
    span = 10;
end

hn = span + intercept;
mn = intercept;
he = rho_a*span;
me = rho_a*rho_b*span/(1+rho_b);
lossmat = loss_matrix(he,hn,me,mn);

end