%% Likelihood function calculated via Durbin Levinson
function [l] = lhoodAR2(paramvector,y)
% Paramvector is the vector of parameters
c=paramvector(1);
Phi=paramvector(2:3);
sigmasq=paramvector(4);

% y is the timeseries - sample data
n=length(y); 
%calculate ACVF
A=ACVF_AR2(n-1,Phi,sigmasq);
% Starting values for D-L recursion & lhood
v=A(1);
l=-0.5*(log(v)+((y(1)-c)^2)/v);

a=A(2)/A(1);
v=A(1)*(1-a^2);
p=c+a*y(1);
l=l-0.5*(log(v)+((y(2)-p)^2)/v);
% D-L recursion & calculation of lhood
for   i=2:n-1
    a_n=(A(i+1)-a'*A(i:-1:2))/v;
    a_1=a-a_n*[flipud(a)];
    a=[a_1;a_n];
    v=v*(1-(a_n)^2);
    p=c+a'*y(i:-1:1);
    l=l-0.5*(log(v)+((y(i+1)-p)^2)/v);
end

l=-(l-n/(2*pi));
end