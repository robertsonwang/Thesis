%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adv. Econometric Methods III                                          %
%      Phillips-Perron test for case 2                                       %
%  h - result of ADF test. 0: Null, 1: Alternative
%  t - adjusted t-stat
%  c - critical value
%
% This function performs a Phillips-Perron test of whether or not a time
% series has a unit root. Outputs are as above. The inputs are:
% y - time series
% q - highest order of sample autocovariance used in calculating NW estimator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[h,t,c,b,SE]=PPtest(y,q)

% y - time series
% q - # lags for NW estimator

% Generate data matrix

T=length(y);

L=lagmatrix(y,1);

X=ones(T, 2);

X(:,1)=L;

X=X(2:end,:);
y=y(2:end);
k=length(X(1,:));
T=length(y);

%perform regression

b=regress(y,X);

resid = y-X*b;

s_sq = (1/(T-k))*resid'*resid;

V = s_sq*inv(X'*X);

SE = diag(sqrt(V));

% calculate \gamma_0
gamma_zero = (1/T)*resid'*resid;

%recursive calculation of Newey-West estimator.
NW =gamma_zero;

for i = 1:q
  
sum=0;
for j=i+1:T
   sum = resid(j)*resid(j-i);
end

NW = NW+(2)*(1-j/(q+1))*(sum/T);

end

%Calculation of modified t-stat
ols_t = (b(1)-1)/SE(1);
t = sqrt((gamma_zero/NW))*ols_t - (1/2)*((NW-gamma_zero)/sqrt(NW))*(T*SE(1))/sqrt(s_sq);

c= DF_Case2_cValue(T); 

h = 0; %can't reject the null (ADF stat)

if t<c
    h = 1; % Reject null of unit root
end


