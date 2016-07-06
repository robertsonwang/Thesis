%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adv. Econometric Methods III                                          %
%      Dickey-Fuller test                                               %
%  h - result of ADF test. 0: Null, 1: Alternative
%  t - t-stat
%  c - critical value
%  F - F stat for joint sig of all lags
%  c_value - critical value for the F test
%   Joint_Sig -  0: null (insignificant)
%
% This function perfoms a Dickey-Fuller/ADF test of whether or not a time
% series has a unit root. Outputs are as above. The inputs are:
%
% y - time series
% lags - number of lagged differences used; 0 for DF
% model - Case 1, 2 or 3, in the terminology used in Hamilton.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function[h,t,c,F,c_value,Joint_Sig,b]= DFTest(y,lags,model)


T=length(y);
L = lagmatrix(y,[0:lags+1]);


if model == 1
 
% Case 1 according to Hamilton. See pg. 528 -  critical value at 5% is
% always -1.95. 


%Generate data matrix
X=zeros(T,lags+1);
X(:,1) = L(:,2);
for j=1:lags
    X(:,j+1)=L(:,j+1)-L(:,j+2);
end
k = length(X(1,:));
start_row=find(isnan(X(:,end)),1,'last')+1; %drop rows with NaN
X=X(start_row:end,:);
y=y(start_row:end);
T=length(y);
%OLS
b=regress(y,X); 

resid = y-X*b;

s_sq = (1/(T-k))*resid'*resid;

V = s_sq*inv(X'*X);

SE = diag(sqrt(V));

% t - statt for coefficient on AR(1) coefficient
t=(b(1)-1)/SE(1);

% F - test for joint significance of all lags. 
if lags > 0 

I_k = eye(k);
R=I_k(2:end,:); 

F = ((R*b)'*inv(R*inv(X'*X)*R')*(R*b))/((k-1)*s_sq);
c_value = finv(0.95,k-1,T-k);

Joint_Sig = 0; %null of joint insignificance

if abs(F)>c_value
    Joint_Sig=1; %H_1 for F_test
end

else
   
 F = (b(1)-0)/SE(1); %t test is same as F test for one restriction
    c_value = tinv(0.975,T-k);
   
    Joint_Sig = 0; %null of joint insignificance

if abs(F)>c_value
    Joint_Sig=1; %H_1
end
end


c=-1.95;
h = 0; %can't reject the null (ADF stat)
if t<c
    h = 1; %can reject the null of a unit root
end

elseif model == 2
    
    % This is for case 2 in Hamilton's terminology.
    
%generate data matrix
X=ones(T,lags+2);
X(:,1) = L(:,2);
for j=1:lags
    X(:,j+1)=L(:,j+1)-L(:,j+2);
end
k = length(X(1,:));
start_row=find(isnan(X(:,end-1)),1,'last')+1; % drop NaN rows
X=X(start_row:end,:);
y=y(start_row:end);
T=length(y);
%OLS
b=regress(y,X); 

resid = y-X*b;

s_sq = (1/(T-k))*resid'*resid;

V = s_sq*inv(X'*X);

SE = diag(sqrt(V));

% t - stat for coefficient on AR(1) coefficient
t=(b(1)-1)/SE(1);

% F - test for joint significance of all lags. 
if lags > 0 

I_k = eye(k);
R=I_k(2:end-1,:); 

F = ((R*b)'*inv(R*inv(X'*X)*R')*(R*b))/((k-1)*s_sq);
c_value = finv(0.95,k-1,T-k);

Joint_Sig = 0; %null of joint insignificance

if abs(F)>c_value
    Joint_Sig=1; %H_1
end

else %If lags  = 0
       
 F = (b(1)-0)/SE(1); %t test is same as F test for one restriction
    c_value = tinv(0.975,T-k);
   
    Joint_Sig = 0; %null of joint insignificance

if abs(F)>c_value
    Joint_Sig=1; %H_1
end
end

c= DF_Case2_cValue(T); 

h = 0; %can't reject the null (ADF stat)

if t<c
    h = 1; %Can reject the null of a unit root
end
   
    

elseif model == 3  % The following is for case 3, in Hamilton's definition. 
    % that is, the true process is y_t = \alpha + y_{t-1} +u_t, with
    % \alpha \neq 0.
    % Hamilton p. 487 gives that in this case, the usual t and F tests
    % have the usual t and F distributions.


%generate data matrix
X=ones(T,lags+2);
X(:,1) = L(:,2);
for j=1:lags
    X(:,j+1)=L(:,j+1)-L(:,j+2);
end
k = length(X(1,:));
start_row=find(isnan(X(:,end-1)),1,'last')+1; %drop NaN rows
X=X(start_row:end,:);
y=y(start_row:end);
T=length(y);
%OLS
b=regress(y,X); 

resid = y-X*b;

s_sq = (1/(T-k))*resid'*resid;

V = s_sq*inv(X'*X);

SE = diag(sqrt(V));

% t - stat for coefficient on AR(1) coefficient
t=(b(1)-1)/SE(1);

% F - test for joint significance of all lags. 
if lags > 0 

I_k = eye(k);
R=I_k(2:end-1,:); 

F = ((R*b)'*inv(R*inv(X'*X)*R')*(R*b))/((k-1)*s_sq);
c_value = finv(0.95,k-1,T-k);

Joint_Sig = 0; %null of joint insignificance

if abs(F)>c_value
    Joint_Sig=1; % H_1
end

else
       
 F = (b(1)-0)/SE(1); %t test is same as F test for one restriction
    c_value = tinv(0.975,T-k);
   
    Joint_Sig = 0; %null of joint insignificance

if abs(F)>c_value
    Joint_Sig=1; % H_1
end
end

c= tinv(0.975, T-k); 

h = 0; %can't reject the null (ADF stat)

if abs(t)>c
    h = 1; % Can reject the null of a unit root.
end

end






