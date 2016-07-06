%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adv. Econometric Methods III                                          %
%      Empirical Homework 2                                             %
%                                                                       %
%                                                                       %
% Team 3:                                                               %
% Suleman Dawood, Bjarni Einarsson, Adam Lee & Robertson Wang           %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping
clear all
close all

clc
%% Data prep - Quarterly

[data0, names] = xlsread('AllData.xlsx'); % load data
[data1, names] = xlsread('AllDataInflation.xls');
ur = data0(:,5);    % UK unemployment rate
inf = data1(:,5);    % UK inflation rate
year = data0(:,1);
month = data0(:,2);


T = length(ur);    % Number of quarters 
%TT_inf=(month_inf(1)/12)+year_inf(1):(1/12):(month_inf(end)/12)+year_inf(end);
%
TT=(month(1)/12)+year(1):(1/12):(month(end)/12)+year(end);


%% Data Prep - Quarterly

k=zeros(T,1);
for i=1:T
    if (month(i) == 2)|(month(i) == 5)|(month(i) == 8)|(month(i) == 11);
        k(i)=1;
end
end

ur_q = ur(k==1);
inf_q = inf(k==1);
year_q=year(k==1);
month_q=month(k==1);
T_q=length(ur_q);
TT_q=(month_q(1)/12)+year_q(1):(3/12):(month_q(end)/12)+year_q(end);


%% Data Exploration - Monthly
% 1. Plot of the time series, inflation vs. unemployment
figure; hold on
a1 = plot(TT,[ur inf]); 
legend('Unemployment','Inflation')

% 2a. Plot of the autocorrelation function for inflation 
figure
autocorr(inf)
% 2b. Plot of the autocorrelation function for unemployment 
figure
autocorr(ur)

% 3a. Plot of the periodogram for inflation 
figure
periodogram(inf)
% 3b. Plot of the periodogram function for unemployment 
figure
periodogram(ur)


%% Data Exploration - Quarterly
% 5. Plot of the time series, inflation vs. unemployment
figure; hold on
a1 = plot(TT_q,[ur_q inf_q]); 
legend('Unemployment','Inflation')

% 2a. Plot of the autocorrelation function for inflation 
figure
autocorr(inf_q)
% 2b. Plot of the autocorrelation function for unemployment 
figure
autocorr(ur_q)

% 3a. Plot of the periodogram for inflation 
figure
periodogram(inf_q)
% 3b. Plot of the periodogram function for unemployment 
figure
periodogram(ur_q)

%% 2a - Unit Root Testing, DF

% a) Dickey-Fuller for Inflation, Monthly

J=0;
lags=1;
while J==0

[h_inf_m,t_inf_m,c_inf_m,f_inf_m,c2_inf_m,J,b_inf_m] = DFTest(inf,lags,2);

lags=lags+1;

end
lags_inf_m = lags-1;

% b) Dickey-Fuller for UR, Monthly

J=0;
lags=1;
while J==0

[h_ur_m,t_ur_m,c_ur_m,f_ur_m,c2_ur_m,J,b_ur_m] = DFTest(ur,lags,2);

lags=lags+1;

end
lags_ur_m = lags-1;

% c) Dickey-Fuller for Inflation, Quarterly

J=0;
lags=1;
while J==0

[h_inf_q,t_inf_q,c_inf_q,f_inf_q,c2_inf_q,J,b_inf_q] = DFTest(inf_q,lags,2);

lags=lags+1;

end
lags_inf_q = lags-1;

% d) Dickey-Fuller for UR, Quarterly

J=0;
lags=1;
while J==0

[h_ur_q,t_ur_q,c_ur_q,f_ur_q,c2_ur_q,J,b_ur_q] = DFTest(ur_q,lags,2);

lags=lags+1;

end
lags_ur_q = lags-1;

%% 2b - Unit Root Testing, PP

h_m=4; %Lag length for NW estimator, monthly
h_q=3; %Lag length for NW estimator, quarterly

% a) Phillips - Perron for Inflation, Monthly
[h_inf_m_pp,t_inf_m_pp,c_inf_m_pp,b_inf_m_pp,SE_inf_m_pp] = PPtest(inf,h_m);
% b) Phillips - Perron for Unemployment, Monthly
[h_ur_m_pp,t_ur_m_pp,c_ur_m_pp,b_ur_m_pp,SE_ur_m_pp] = PPtest(ur,h_m);
% c) Phillips - Perron for Inflation, Quarterly
[h_inf_q_pp,t_inf_q_pp,c_inf_q_pp,b_inf_q_pp,SE_inf_q_pp] = PPtest(inf_q,h_q);
% d) Phillips - Perron for Unemployment, Quarterly
[h_ur_q_pp,t_ur_q_pp,c_ur_q_pp,b_ur_q_pp,SE_ur_q_pp] = PPtest(ur_q,h_q);

%% 3 - Descriptive VAR analysis
% - Lag selection using BIC
L = 4;              % Lag length
BIC = NaN(L,1);
Y = [ur_q inf_q];
n = size(Y,2);
for j = 1:L
    X = [ones(size(Y,1),1)];
    for i = 1:j
        X= [X lagger(Y,i)];
    end
    Y = Y(j+1:end,:);
    X = X(j+1:end,:);
    Ts = size(Y,1);
    
    bols = inv(X'*X)*(X'*Y);
    Sigma = ((Y-X*bols)'*(Y-X*bols))/Ts;
    BIC(j,1) = log(det(Sigma)) + log(Ts)/Ts*j*n^2;
end

% - Estimate preferred model
minbic = min(BIC);
BIC_min = BIC==minbic;
L=find(BIC_min);
Y = [ur_q inf_q];
X = [ones(size(Y,1),1)];
for i = 1:L
    X= [X lagger(Y,i)];
end
Y = Y(L+1:end,:);
X = X(L+1:end,:);
Ts = size(Y,1);
k=size(X,2);
bols = inv(X'*X)*(X'*Y);
Sigma = ((Y-X*bols)'*(Y-X*bols))/Ts;


%% Granger Causality


% Equation for ur
R1 = [0 0 1 0 0;0 0 0 0 1];
F_ur = (R1*bols(:,1))'*inv(R1*inv(X'*X)*R1')*(R1*bols(:,1))/(Sigma(1,1)); %#R is 1


% Equation for inf
R2 = [0 1 0 0 0;0 0 0 1 0];
F_inf = (R2*bols(:,2))'*inv(R2*inv(X'*X)*R2')*(R2*bols(:,2))/(Sigma(2,2)); %#R is 1

% Critical Value

Granger_cValue = finv(0.95,1,Ts-k);

if F_ur>Granger_cValue
    Granger_F_ur = 1; %Alternative
else Granger_F_ur=0; %Null - inf does not Granger-Cause ur
end

if F_inf>Granger_cValue
    Granger_F_inf = 1; %Alternative
else Granger_F_inf=0; %Null - ur does not Granger-Cause inf
end
    


%% Impulse response functions
% We note that the preferred model has p = 2.
% Transform AR(2) into AR(1) companion form

Phi1 = [bols(2,1),bols(3,1); bols(2,2), bols(3,2)];
Phi2 = [bols(4,1),bols(5,1); bols(4,2), bols(5,2)];

A = [Phi1, Phi2; eye(2,2), zeros(2,2)];
mu = [bols(1,1);bols(1,2)];
Q = chol(Sigma, 'lower');

C = [Q ; zeros(2,2)];


h=100;

   % h period IRF based on a shock to ur
   shock_ur = [1;0];
   C1 = C*shock_ur;
   
   
   IRF_ur  = zeros(4,h);
   
   for i=1:h
       IRF_ur(:,(i-1)+1) = A^(i-1)*C*shock_ur;
   end
   ur_irf_urshock = IRF_ur(1,:);  
   inf_irf_urshock = IRF_ur(2,:); 
   
   figure
   subplot(2,1,1)
   plot(ur_irf_urshock)
   title('Impulse response of unemployment to an unemployment shock')
   subplot(2,1,2)
   plot(inf_irf_urshock)
   title('Impulse response of inflation to an unemployment shock')
   %saveas(gcf,'irf_urshock.png')
   
   % h period IRF based on a shock to inf
   shock_inf = [0;1];
   C1 = C*shock_inf;
   
   IRF_inf  = zeros(4,h);
   
   for i=1:h
       IRF_inf(:,(i-1)+1) = A^(i-1)*C*shock_inf;
   end
   ur_irf_infshock = IRF_inf(1,:);  
   inf_irf_infshock = IRF_inf(2,:);  
   
   figure
   subplot(2,1,1)
   plot(ur_irf_infshock)
   title('Impulse response of unemployment to an inflation shock')
   subplot(2,1,2)
   plot(inf_irf_infshock)
   title('Impulse response of inflation to an inflation shock')
   %saveas(gcf,'irf_infshock.png')


   
   %% Variance Decomposition
   
   VD_ur=zeros(h,2);
   for i=1:h
       VD_ur(i,1)=FEVD(i,1,1,Sigma,A);
       VD_ur(i,2)=FEVD(i,1,2,Sigma,A);
   end
   
   VD_inf=zeros(h,2);
   for i=1:h
       VD_inf(i,1)=FEVD(i,2,1,Sigma,A);
       VD_inf(i,2)=FEVD(i,2,2,Sigma,A);
   end
   
   figure 
   bar(VD_ur,'stacked')
   axis('tight')
   legend('FEV due to shock to unemployment ', 'FEV due to shock to inflation')
   title('Forecast error variance decomposition for Unemployment')
   %saveas(gcf,'FEVD_ur.png')
   figure
   bar(VD_inf,'stacked')
   axis('tight')
   legend('FEV due to shock to unemployment', 'FEV due to shock to inflation')
   title('Forecast error variance decomposition for Inflation')
   %saveas(gcf,'FEVD_inf.png')


%% 4 - Time Varying Paramters
    %% Estimation
options = optimset('PlotFcns',@optimplotx,'DISPLAY','iter');
AA = [0,-1,0;0,0,-1];
BB= [-0.000000000001;-0.000000000001];
[param,~,~,~,~,~,hess] = fmincon(@(param) kfloglik(param,[ur inf]), [0.1 0.1 0.1],AA,BB,[],[],-Inf,Inf,[],options);
separam = sqrt(diag(inv(hess)));


    %% Kalman output
[ap, pp, af, pf, as, ps,r, N] = kalmansmooth(param,[ur inf]);


    %% Figures
figure
plot([af, af+1.96*sqrt(pf), af-1.96*sqrt(pf)])

figure 
plot([as(2:end), as(2:end)+1.96*real(sqrt(ps(2:end))), as(2:end)-1.96*real(sqrt(ps(2:end)))])


