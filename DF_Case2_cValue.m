%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adv. Econometric Methods III                                          %
%      Critical Values for Dickey-Fuller test  in case 2                %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function c = DF_Case2_cValue(T)

if T <= 25
    c = -3;
elseif T<=50
    c = -2.93;
elseif T<=100
    c=-2.89;
elseif T<=250
    c=-2.88;
elseif T<=500
    c=-2.87;
else c=-2.86;
end
     