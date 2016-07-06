%% Implements the Asymptotic test in Sec 1.1. of Diebold Mariano (1995)
%
% INPUTS
%   FE1: Forecast error 1
%   FE2: Forecast error 2
%   steps: how many steps forward the forecast is
%
% OUTPUTS
%   pval: the p-value of the test of the null of no difference
%
function [DM, pval] = DieboldMariano(FE1,FE2,steps)

%% 
T = size(FE1,1);  % number of observations
d = FE1-FE2;    % forecast error difference
dbar = mean(d);
gamma0 = var(d);
if steps > 1    % correction for autocovariance if more than 1 step forecast
    gamma = zeros(steps-1,1);
    for i = 1:steps-1
        covbar = cov(d(1+i:T),d(1:T-i));
        gamma(i) = covbar(2);
    end
    dvar = gamma0 + 2*sum(gamma);
else
    dvar = gamma0;
end
DM = dbar /sqrt((1/T)*dvar);        % Test statistic
pval = 2*(1-cdf('Normal',abs(DM),0,1));  % p-value


    