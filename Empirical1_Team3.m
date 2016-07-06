%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adv. Econometric Methods III                                          %
%      Empirical Homework 1                                             %
%        - The objective of this problem set is forecasting             %
%          unemployment in the United Kingdom using ARMA models         %
%                                                                       %
% Team 3:                                                               %
% Suleman Dawood, Bjarni Einarsson, Adam Lee & Robertson Wang           %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping
clear all
close all


%% Data prep
[data0, names] = xlsread('AllData.xlsx'); % load data
ur = data0(:,5);    % UK unemployment rate
year = data0(:,1);
month = data0(:,2);
T = size(data0,1);    % Number of periods in sample
TT = 1986.3333:(1/12):2016;
% Crisis start and end dates
startdate = [1973.25, 1974.5, 1979.25, 1990.25, 2008]';
enddate = [1974, 1975.5, 1981, 1991.5, 2009.25]';


%% Data Exploration
% 1. Plot of the time series
plot(TT',ur)
curax = axis;
title('UK unemployment rate')

% 2. Plot with recession indicator
figure
hold on
a4 = area([startdate(4) enddate(4)], [curax(4) curax(4)],curax(3));
a4.FaceColor = [0.8,0.8,0.8];
a4.EdgeColor = [0.8,0.8,0.8];
a5 = area([startdate(5) enddate(5)], [curax(4) curax(4)],curax(3));
a5.FaceColor = [0.8,0.8,0.8];
a5.EdgeColor = [0.8,0.8,0.8];
plot(TT',ur)
hold off
title('UK unemployment rate and recessions')

% 3. See .pdf

% 4. Autocorrelation function
figure
autocorr(ur)

% 5. Partial autocorrelation function
figure
parcorr(ur)

% 6. Histogram
figure
hist(ur)
title('Histogram of UK unemplyment rate')

%% 2 Models and Estimation

% 2.1 ARMA models
%    AR(2)
startvalAR2 = [0.1,0.1,0.1,0.1];
Aar2=[0,1,1,0;0,-1,1,0;0,0,1,0;0,0,-1,0];
bar2=[0.99999;0.99999;0.99999;0.99999];
[thetaAR2,fvalAR2,~,~,~,~,hessian] = fmincon(@(theta)lhoodAR2(theta,ur),startvalAR2,Aar2,bar2);
errAR2 =sqrt(diag(inv(hessian)));
resAR2 = zeros(T,1);
for i = 3:T
    resAR2(i) = ur(i)-thetaAR2(1) - thetaAR2(2)*ur(i-1) - thetaAR2(3)*ur(i-2);
end
figure
autocorr(resAR2)
figure
parcorr(resAR2)

%    ARMA(1,1)
startvalARMA11 = [0.1,0.1,0.1,0.1];
Aarma11=[0,1,0,0;0,-1,0,0];
barma11=[1;1];
[thetaARMA11,fvalARMA11,~,~,~,~,hessian] = fmincon(@(theta)lhoodARMA11(theta,ur),startvalARMA11,Aarma11,barma11);
errARMA11=sqrt(diag(inv(hessian)));
resARMA11 = zeros(T,1);
for i = 2:T
    resARMA11(i) = ur(i)-thetaARMA11(1)-thetaARMA11(2)*ur(i-1)-thetaARMA11(3)*resARMA11(i-1);
end
figure
autocorr(resARMA11)
figure
parcorr(resARMA11)


% 2.2 Model selection
maxp = 2;
maxq = 2;
AIC = NaN(maxp+1,maxq+1);
BIC = NaN(maxp+1,maxq+1);
AICC = NaN(maxp+1,maxq+1);
loglik = NaN(maxp+1,maxq+1);
loglik(3,1) = fvalAR2;
loglik(2,2) = fvalARMA11;
% AR(1)
A=[0,1,0;0,-1,0];
b=[0.99999;0.99999];
theta0=[0.1,0.1,0.1];
[~,fvalAR1,~,~,~,~] = fmincon(@(theta)lhoodAR1(theta,ur),theta0,A,b);
loglik(2,1) = fvalAR1;
% MA(1)
theta0=[0.1,0.1,0.1];
[~,fvalMA1,~,~,~,~] = fminunc(@(theta)lhoodMA1(theta,ur),theta0);
loglik(1,2) = fvalMA1;
% MA(2)
theta0=[0.1,0.1,0.1,0.1];
[~,fvalMA2,~,~,~,~] = fminunc(@(theta)lhoodMA2(theta,ur),theta0);
loglik(1,3) = fvalMA2;
% ARMA(1,2)
theta0=[0.1,0.1,0.1,0.1,0.1];
A=[0,1,0,0,0;0,-1,0,0,0];
b=[1;1]; 
[~,fvalARMA12,~,~,~,~,~] = fmincon(@(theta)lhoodARMA12(theta,ur),theta0,A,b);
loglik(2,3) = fvalARMA12;
% ARMA(2,1)
A=[0,1,1,0,0;0,-1,1,0,0;0,0,1,0,0;0,0,-1,0,0];
b=[0.99999;0.99999;0.99999;0.99999];
theta0=[0.1,0.1,0.1,0.1,0.1];
[~,fvalARMA21,~,~,~,~,~] = fmincon(@(theta)lhoodARMA21(theta,ur),theta0,A,b);
loglik(3,2) = fvalARMA21;
% ARMA(2,2)
for p = 0:maxp
    for q = 0:maxq
        AIC(p+1,q+1) = 2*loglik(p+1,q+1) + 2*(p+q+1);
        BIC(p+1,q+1) = 2*loglik(p+1,q+1) + (p+q+1)*(log(T)/T);
        AICC(p+1,q+1) = 2*loglik(p+1,q+1) + 2*T*(p+q+1)/(T-p-q-2);
    end
end
% All criteria agree on AR(2) except BIC (ARMA(2,1)), AR(2) chosen for
% parsimony

%% 3 Forecasting
hh = 6; % longest forecast horizon
strtpt = 226; % Jan 2005 = obs 226
endpt = T;
yhatAR2 = zeros(endpt-strtpt,2);
yhatARMA = zeros(endpt-strtpt,2);
wait = waitbar(0,'Forecasting algorithm is running...');
for i = 0:endpt-strtpt
    % Adjust sample
    smpl = ur(1:strtpt+i); 
    % Estimate models
    [thetaAR2,~,~,~,~,~,~] = fmincon(@(theta)lhoodAR2(theta,smpl),startvalAR2,Aar2,bar2);
    [thetaARMA11,~,~,~,~,~,~] = fmincon(@(theta)lhoodARMA11(theta,smpl),startvalARMA11,Aarma11,barma11);
    resid = zeros(size(smpl,1),1);
    for j = 2:size(smpl,1)
        resid(j) = ur(j)-thetaARMA11(1)-thetaARMA11(2)*ur(j-1)-thetaARMA11(3)*resid(j-1);
    end
    % Calculate forecasts
    %   AR(2)
    temp = zeros(hh+2,1);
    temp(1:2) = smpl(end-1:end);
    for j = 3:length(temp)
        temp(j) = thetaAR2(1:end-1)*[1;flipud(temp(j-2:j-1))];
    end
    yhatAR2(i+1,:) = [temp(3) temp(end)];
    
    %   ARMA(1,1)
    temp = zeros(hh+1,2);
    temp(1,1) = smpl(end);
    temp(1,2) = resid(end);
    for j = 2:length(temp)
        temp(j) = thetaARMA11(1:end-1)*[1;temp(j-1,1);temp(j-1,2)];
    end
    yhatARMA(i+1,:) = [temp(2) temp(end)];
    
    % Waitbar progress
    if mod(i,5)==0
        waitbar(i/(endpt-strtpt), wait);
    end

end

% 1. Plot the forecasts
figure
plot(TT',ur, TT',[ur(1:strtpt);yhatAR2(1:end-1,1)],'--', TT',[ur(1:strtpt+hh-1);yhatAR2(1:end-6,2)],':')
title('AR(2)');
 
figure
plot(TT',ur, TT',[ur(1:strtpt);yhatARMA(1:end-1,1)],'--', TT',[ur(1:strtpt+hh-1);yhatARMA(1:end-6,2)],':')
title('ARMA(1,1)')

% 2. MSFE & MAFE
% AR(2)
FE1sAR = ur(strtpt+1:end)-yhatAR2(1:end-1,1); 
FE6sAR = ur(strtpt+hh:end)-yhatAR2(1:end-6,2);
MSFE1sAR = mean(FE1sAR.^2);
MSFE6sAR = mean(FE6sAR.^2);
MAFE1sAR = mean(abs(FE1sAR));
MAFE6sAR = mean(abs(FE6sAR));
 
% ARMA(1,1)
FE1sARMA = ur(strtpt+1:end)-yhatARMA(1:end-1,1); 
FE6sARMA = ur(strtpt+hh:end)-yhatARMA(1:end-6,2);
MSFE1sARMA = mean(FE1sARMA.^2);
MSFE6sARMA = mean(FE6sARMA.^2);
MAFE1sARMA = mean(abs(FE1sARMA));
MAFE6sARMA = mean(abs(FE6sARMA));

% 3. See .pdf


%% 4 Forecast evaluation
% 
[DMsq, pvalsq] = DieboldMariano(FE6sAR.^2,FE6sARMA.^2,6);
[DMabs, pvalabs] = DieboldMariano(abs(FE6sAR),abs(FE6sARMA),6);




