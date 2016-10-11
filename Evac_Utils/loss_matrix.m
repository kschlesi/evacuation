function lossmat = loss_matrix(he,hn,me,mn)
% returns the loss matrix function handle with the four constants included.
% the returned function handle takes two boolean inputs, didHit and didEvac
% INPUTS: losses on (hit/evac, hit/noevac, miss/evac, miss/noevac)

% defaults as used in experiment
if ~nargin
    hn = 10;
    he = 6;
    me = 2;
    mn = 0;
end

lossmat = @(didHit,didEvac) hn*(didHit.*~didEvac) + ...
                            he*(didHit.*didEvac) + ...
                            mn*(~didHit.*~didEvac) + ...
                            me*(~didHit.*didEvac) ;                     
end