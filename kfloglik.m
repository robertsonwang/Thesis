%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adv. Econometric Methods III                                          %
%      Empirical Homework 2                                             %
%                                                                       %
% Team 3:                                                               %
% Suleman Dawood, Bjarni Einarsson, Adam Lee & Robertson Wang           %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ll = kfloglik(param,data)

phi = param(1);
H = param(2);
Q = param(3);

y = data(2:end,1);
n = size(y,1);
Z = lagger(data,1);
yl = Z(2:end,1);
Z = Z(2:end,2);
% Initialisation
ll = 0;
s = zeros(n,1);
s2 = zeros(n,1);
a10 = zeros(1,1);
P10 = eye(1)*10^7;
% Matrices
T = eye(1);
R = [1];

for i = 1:n
    v   = y(i,:)-phi*yl(i,:)-Z(i,:)*a10;
    F   = Z(i,:)*P10*Z(i,:)' + H; 
    a11 = a10 + P10*Z(i,:)'*inv(F)*v;
    P11 = P10 - P10*Z(i,:)'*inv(F)*Z(i,:)*P10;
    a10 = T*a11;
    P10 = T*P10*(T-T*P10*Z(i,:)'*inv(F)*Z(i,:))'+R*Q*R';
    % Calculate negative log likelihood
    if i > 1
        ll = ll + 0.5*(log(2*pi) + log(F)+v^2/F);
    end
    
    s(i,:) = a11;
    s2(i,:,:) = P11;
end
    
