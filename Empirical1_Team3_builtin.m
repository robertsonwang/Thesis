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
AR2 = arima(2,0,0);
[estAR2, covAR2, llAR2,~] = estimate(AR2,ur);
[resAR2,~] = infer(estAR2,ur);
figure
autocorr(resAR2)
figure
parcorr(resAR2)

%    ARMA(1,1)
ARMA11 = arima(1,0,1);
[estARMA11,covARMA11, llARMA11,~] = estimate(ARMA11,ur);
[resARMA11,~] = infer(estARMA11,ur);
figure
autocorr(resARMA11)
figure
parcorr(resARMA11)


% 2.2 Model selection
maxp = 5;
maxq = 5;
AIC = NaN(maxp+1,maxq+1);
BIC = NaN(maxp+1,maxq+1);
AICC = NaN(maxp+1,maxq+1);
loglik = NaN(maxp+1,maxq+1);
for q = 1:maxq
    tempmdl = arima(0,0,q);
    [~,~,llik,~] = estimate(tempmdl,ur);
    AIC(1,q+1) = -2*llik + 2*(q+1);
    BIC(1,q+1) = -2*llik + (q+1)*(log(T)/T);
    AICC(1,q+1) = -2*llik + 2*T*(q+1)/(T-q-2);
end
for p = 1:maxp
    for q = 0:maxq
        tempmdl = arima(p,0,q);
        [~,~,llik,~] = estimate(tempmdl,ur);
        AIC(p+1,q+1) = -2*llik + 2*(p+q+1);
        BIC(p+1,q+1) = -2*llik + (p+q+1)*(log(T)/T);
        AICC(p+1,q+1) = -2*llik + 2*T*(p+q+1)/(T-p-q-2);
    end
end

% Preferred model - ARMA(2,1)
ARMA21 = arima(2,0,1);
[estARMA21,covARMA21, llARMA21,~] = estimate(ARMA21,ur);


%% 3 Forecasting
hh = 6; % longest forecast horizon
strtpt = 226; % Jan 2005 = obs 226
endpt = T;
yhatAR2 = zeros(endpt-strtpt,2);
yhatARMA11 = zeros(endpt-strtpt,2);
yhatARMA21 = zeros(endpt-strtpt,2);
wait = waitbar(0,'Forecasting algorithm is running...');
for i = 0:endpt-strtpt
    % Adjust sample
    smpl = ur(1:strtpt+i); 
    % Estimate models
    %  AR(2)
    mdlAR2 =arima(2,0,0);
    estmdlAR2 = estimate(mdlAR2,smpl);
       
    %  ARMA(1,1)
    mdlARMA11 = arima(1,0,1);
    estmdlARMA11 = estimate(mdlARMA11,smpl);
    [resARMA11,~]=infer(estmdlARMA11,smpl);
    
    %  ARMA(2,1)
    mdlARMA21 = arima(2,0,1);
    estmdlARMA21 = estimate(mdlARMA21,smpl);
    [resARMA21,~] = infer(estmdlARMA21,smpl);
    
    % Calculate forecasts
    %   AR(2)
    temp = zeros(hh+2,1);
    temp(1:2) = smpl(end-1:end);
    for j = 3:length(temp)
        temp(j) = [estmdlAR2.Constant cell2mat(estmdlAR2.AR)]*[1;flipud(temp(j-2:j-1))];
    end
    yhatAR2(i+1,:) = [temp(3) temp(end)];
    
    %   ARMA(1,1)
    temp = zeros(hh+1,2);
    temp(1,1) = smpl(end);
    temp(1,2) = resARMA11(end);
    for j = 2:length(temp)
        temp(j) = [estmdlARMA11.Constant cell2mat(estmdlARMA11.AR) cell2mat(estmdlARMA11.MA)]*[1;temp(j-1,1);temp(j-1,2)];
    end
    yhatARMA11(i+1,:) = [temp(2,1) temp(end,1)];
    
    %   ARMA(2,1)
    temp = zeros(hh+2,2);
    temp(1:2,1) = smpl(end-1:end);
    temp(2,2) = resARMA21(end);
    for j = 3:length(temp)
        temp(j) = [estmdlARMA21.Constant cell2mat(estmdlARMA21.AR) cell2mat(estmdlARMA21.MA)]*[1;flipud(temp(j-2:j-1,1));temp(j-1,2)];
    end
    yhatARMA21(i+1,:) = [temp(3,1) temp(end,1)];
    
    % Waitbar progress
    if mod(i,5)==0
        waitbar(i/(endpt-strtpt), wait);
    end

end
close(wait);

% 1. Plot the forecasts
figure
plot(TT',ur, TT',[ur(1:strtpt);yhatAR2(1:end-1,1)],'--', TT',[ur(1:strtpt+hh-1);yhatAR2(1:end-6,2)],':')
title('AR(2)');
 
figure
plot(TT',ur, TT',[ur(1:strtpt);yhatARMA11(1:end-1,1)],'--', TT',[ur(1:strtpt+hh-1);yhatARMA11(1:end-6,2)],':')
title('ARMA(1,1)')

figure
plot(TT',ur, TT',[ur(1:strtpt);yhatARMA21(1:end-1,1)],'--', TT',[ur(1:strtpt+hh-1);yhatARMA21(1:end-6,2)],':')
title('ARMA(2,1)')

% 2. MSFE & MAFE
% AR(2)
FE1sAR = ur(strtpt+1:end)-yhatAR2(1:end-1,1); 
FE6sAR = ur(strtpt+hh:end)-yhatAR2(1:end-6,2);
MSFE1sAR = mean(FE1sAR.^2);
MSFE6sAR = mean(FE6sAR.^2);
MAFE1sAR = mean(abs(FE1sAR));
MAFE6sAR = mean(abs(FE6sAR));
 
% ARMA(1,1)
FE1sARMA11 = ur(strtpt+1:end)-yhatARMA11(1:end-1,1); 
FE6sARMA11 = ur(strtpt+hh:end)-yhatARMA11(1:end-6,2);
MSFE1sARMA11 = mean(FE1sARMA11.^2);
MSFE6sARMA11 = mean(FE6sARMA11.^2);
MAFE1sARMA11 = mean(abs(FE1sARMA11));
MAFE6sARMA11 = mean(abs(FE6sARMA11));

% ARMA(2,1)
FE1sARMA21 = ur(strtpt+1:end)-yhatARMA21(1:end-1,1); 
FE6sARMA21 = ur(strtpt+hh:end)-yhatARMA21(1:end-6,2);
MSFE1sARMA21 = mean(FE1sARMA21.^2);
MSFE6sARMA21 = mean(FE6sARMA21.^2);
MAFE1sARMA21 = mean(abs(FE1sARMA21));
MAFE6sARMA21 = mean(abs(FE6sARMA21));

% 3. See .pdf


%% 4 Forecast evaluation
% AR2 vs ARMA11
[DMsq1, pvalsq1] = DieboldMariano(FE6sAR.^2,FE6sARMA11.^2,6);
[DMabs1, pvalabs1] = DieboldMariano(abs(FE6sAR),abs(FE6sARMA11),6);

% AR2 vs ARMA21
[DMsq2, pvalsq2] = DieboldMariano(FE6sAR.^2,FE6sARMA21.^2,6);
[DMabs2, pvalabs2] = DieboldMariano(abs(FE6sAR),abs(FE6sARMA21),6);

% ARMA11 vs ARMA21
[DMsq3, pvalsq3] = DieboldMariano(FE6sARMA11.^2,FE6sARMA21.^2,6);
[DMabs3, pvalabs3] = DieboldMariano(abs(FE6sARMA11),abs(FE6sARMA21),6);


