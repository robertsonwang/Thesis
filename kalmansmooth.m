%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adv. Econometric Methods III                                          %
%      Empirical Homework 2                                             %
%                                                                       %
% Team 3:                                                               %
% Suleman Dawood, Bjarni Einarsson, Adam Lee & Robertson Wang           %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ap, pp, af, pf, as, ps,r,N] = kalmansmooth(param,data)
phi = param(1);
H = param(2);
Q = param(3);

y = data(2:end,1);
n = size(y,1);
Z = lagger(data,1);
yl = Z(2:end,1);
Z = Z(2:end,2);

ap = zeros(n,1);
pp = zeros(n,1);
af = zeros(n,1);
pf = zeros(n,1);
as = zeros(n,1);
ps = zeros(n,1);
vv = zeros(n,1);

%% Run Kalman filter
a10 = zeros(1,1);
P10 = eye(1)*10^7;
% Matrices
T = eye(1);
R = [1];

for i = 1:n
    % Store predicted states
    ap(i,:) = a10;
    pp(i,:) = P10;
    % Kalman filter
    v   = y(i,:)-phi*yl(i,:)-Z(i,:)*a10;
    F   = Z(i,:)*P10*Z(i,:)' + H; 
    a11 = a10 + P10*Z(i,:)'*inv(F)*v;
    P11 = P10 - P10*Z(i,:)'*inv(F)*Z(i,:)*P10;
    a10 = T*a11;
    P10 = T*P10*(T-T*P10*Z(i,:)'*inv(F)*Z(i,:))'+R*Q*R';
    % Store filtered states
    af(i,:) = a11;
    pf(i,:) = P11;
    % store stuff for smoother
    vv(i,:) = v;
end


%% State smoothing recursion
r1 = 0;
N1 = 0;
r = zeros(n,1);
N = zeros(n,1);
for i = n:-1:1
    r0 = Z(i,:)'*inv(F)*vv(i) + (T-T*pp(i,:)*Z(i,:)'*inv(F)*Z(i,:))'*r1;
    N0 = Z(i,:)'*inv(F)*Z(i,:) + (T-T*pp(i,:)*Z(i,:)'*inv(F)*Z(i,:))'*N1*(T-T*pp(i,:)*Z(i,:)'*inv(F)*Z(i,:));
    as(i,:) = ap(i,:) + pp(i,:)*r0;
    ps(i,:) = pp(i,:) - pp(i,:)*N0*pp(i,:);
    r(i,:) = r0;
    N(i,:) = N0;
end






