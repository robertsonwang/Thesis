%% Monte Carlo Simulations of the AR(2) and ARMA(1,1) estimation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This .m file performs various monte carlo simulations of our
% routines for estimating AR(2) and ARMA(1,1) models.


%% AR(2) Test 1

n=100;

resultsAR2=zeros(n,4);

for j=1:n
        
    Mdl = arima('AR',{0.8,0.15},'Constant',0,'Variance',0.04);
    y=simulate(Mdl,200);
    
    A=[0,1,1,0;0,-1,1,0;0,0,1,0;0,0,-1,0];
    b=[0.99999;0.99999;0.99999;0.99999];
    theta0=[0.1,0.1,0.1,0.1];
    theta= fmincon(@(theta)lhoodAR2(theta,y),theta0,A,b);
    
    resultsAR2(j,:)=theta;
    j 
end
mean(resultsAR2)
median(resultsAR2)
plot(resultsAR2)


%% AR(2) Test 2

n=200;

resultsAR2=zeros(n,4);

for j=1:n
        
    Mdl = arima('AR',{0.6,0.2},'Constant',0,'Variance',(0.2)^2);
    y=simulate(Mdl,200);
    
    A=[0,1,1,0;0,-1,1,0;0,0,1,0;0,0,-1,0];
    b=[0.99999;0.99999;0.99999;0.99999];
    theta0=[0.1,0.1,0.1,0.1];
    theta= fmincon(@(theta)lhoodAR2(theta,y),theta0,A,b);
    
    resultsAR2(j,:)=theta;
    j
end
mean(resultsAR2)
median(resultsAR2)
plot(resultsAR2)


%% AR(2) Test 3

n=100;

resultsAR2=zeros(n,4);

for j=1:n
        
    Mdl = arima('AR',{0.95,0.03},'Constant',0,'Variance',(0.2)^2);
    y=simulate(Mdl,200);
    
    A=[0,1,1,0;0,-1,1,0;0,0,1,0;0,0,-1,0];
    b=[0.99999;0.99999;0.99999;0.99999];
    theta0=[0.1,0.1,0.1,0.1];
    theta= fmincon(@(theta)lhoodAR2(theta,y),theta0,A,b);
    
    resultsAR2(j,:)=theta;
    j
end
mean(resultsAR2)
median(resultsAR2)
plot(resultsAR2)


%% AR(2) Test 4

n=100;

resultsAR2=zeros(n,4);

for j=1:n
        
    Mdl = arima('AR',{0.99,0.005},'Constant',0,'Variance',(0.2)^2);
    y=simulate(Mdl,200);
    
    A=[0,1,1,0;0,-1,1,0;0,0,1,0;0,0,-1,0];
    b=[0.99999;0.99999;0.99999;0.99999];
    theta0=[0.1,0.1,0.1,0.1];
    theta= fmincon(@(theta)lhoodAR2(theta,y),theta0,A,b);
    
    resultsAR2(j,:)=theta;
    j
end
mean(resultsAR2)
median(resultsAR2)
plot(resultsAR2)


%% ARMA(1,1) test 1

n=100;
resultsARMA11=zeros(n,4);

for k=1:n
        
    Mdl = arima('AR',{0.6},'D',0,'MA',{0.2},'Constant',0,'Variance',0.04);
    x=simulate(Mdl,200);
    
    A=[0,1,0,0;0,-1,0,0]; 
    b=[0.99999;0.99999];
    theta0=[0.1,0.1,0.1,0.1];
    theta=fmincon(@(theta)lhoodARMA11(theta,x),theta0,A,b);
    
    resultsARMA11(k,:)=theta;
    k
end

mean(resultsARMA11)
median(resultsARMA11)
plot(resultsARMA11)

%% ARMA(1,1) test 2

n=100;
resultsARMA11=zeros(n,4);

for k=1:n
        
    Mdl = arima('AR',{0.8},'D',0,'MA',{0.15},'Constant',0,'Variance',0.04);
    x=simulate(Mdl,200);
    
    A=[0,1,0,0;0,-1,0,0]; 
    b=[0.99999;0.99999];
    theta0=[0.1,0.1,0.1,0.1];
    theta=fmincon(@(theta)lhoodARMA11(theta,x),theta0,A,b);
    
    resultsARMA11(k,:)=theta;
    k
end

mean(resultsARMA11)
median(resultsARMA11)
plot(resultsARMA11)


%% ARMA(1,1) test 3

n=100;
resultsARMA11=zeros(n,4);

for k=1:n
        
    Mdl = arima('AR',{0.95},'D',0,'MA',{0.03},'Constant',0,'Variance',0.04);
    x=simulate(Mdl,200);
    
    A=[0,1,0,0;0,-1,0,0]; 
    b=[0.99999;0.99999];
    theta0=[0.1,0.1,0.1,0.1];
    theta=fmincon(@(theta)lhoodARMA11(theta,x),theta0,A,b);
    
    resultsARMA11(k,:)=theta;
    k
end

mean(resultsARMA11)
median(resultsARMA11)
plot(resultsARMA11)


%% ARMA(1,1) test 4

n=100;
resultsARMA11=zeros(n,4);

for k=1:n
        
    Mdl = arima('AR',{0.99},'D',0,'MA',{0.005},'Constant',0,'Variance',0.04);
    x=simulate(Mdl,200);
    
    A=[0,1,0,0;0,-1,0,0]; 
    b=[0.99999;0.99999];
    theta0=[0.1,0.1,0.1,0.1];
    theta=fmincon(@(theta)lhoodARMA11(theta,x),theta0,A,b);
    
    resultsARMA11(k,:)=theta;
    k
end

mean(resultsARMA11)
median(resultsARMA11)
plot(resultsARMA11)
