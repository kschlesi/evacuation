function [params,ML] = ML_fit_beta(qform,points,H,J,startp,varargin)
% This function finds the best-fit parameters, 'params', for a P-parameter 
% function 'qform' at a set of values 'points'. The fit is performed by
% finding the parameters that maximize the likelihood function under the
% assumption that the data fits a beta distribution, as described in Sean's
% paper. Assumes the ML format given in that paper, in terms of measured
% variables 'H' (total number of opportunities to evacuate) and 'J'(total
% number of evacuations observed
%
% INPUTS: qform,  a handle to a P-parameter function that gives the form of
%                 the decision model, taking 2 inputs: the vector of values
%                 'points', and a P-element vector of parameters.
%         points, a vector of Phit values (or values at which to fit the function)
%         H,      total number of evacuation opportunities at each value in
%                 'points'
%         J,      total evacuations observed at each value in 'points'
%         startp, Px1 vector of starting estimate values for parameters
% OUTPUT: params, Px1 vector of best-fit parameters
%
% NOTE. This restricts all final function output values to be <1 (assuming
% they are probabilities) & uses the values 0:0.1:1 to test this condition

qfit = @(pv_) qform(points,pv_);
MLfit = @(pvec_) -1*sum(((sum(H)-sum(J)).*log(1-qfit(pvec_)) + sum(J).*log(qfit(pvec_))));
options = optimoptions(@fminunc,'MaxFunEvals',10000,varargin{:});
[params1,MLerr] = fminunc(MLfit,startp,options);
A = zeros(length(startp));
b = zeros(length(startp),1);

% constrains all prob. as a function of points to be <=1
q_con = @(pvec_)(qform(0:0.1:1,pvec_)-1); 
[theta_con,MLerr_con] = fmincon(MLfit,startp,A,b,A,b,0,100,@(pvec_)q_add_eq(pvec_,q_con));

% restricts to <1
if MLerr<MLerr_con && all(qform(0:0.1:1,params1)<=1)
    params = params1;
    ML = MLerr;
else
    params = theta_con;
    ML = MLerr_con;
end


%figure(1); hold all; plot(Phit,qform(Phit,params));
%xlabel('Phit'); ylabel('evacuation probability'); title('best-fit hill model, ind. trials, shelterSpace=50');
%legend(['\Lambda = ' num2str(params(1)) ', k = ' num2str(params(2)) ', n = ' num2str(params(3))]);
%legend(['\alpha = ' num2str(params(1)) ', \beta = ' num2str(params(2))]);
%hold off;

end