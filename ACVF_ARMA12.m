function ACVF=ACVF_ARMA12(h,Phi,Theta,sigmasq)
% This computes the ACVF at lag h for a ARMA(1,2) model and returns a
% vector of all lags up to h.

if ~(length(Phi)==1) 
    error('Length of coefficient vector is not 1')
end
if ~(length(Theta)==2)
    error('Length of coefficient vector is not 2')
end


 ACVF_0=(sigmasq/(1-Phi^2))*(1+Theta(1)^2+Theta(1)*Phi+Theta(2)^2+Phi*Theta(1)*Theta(2)+Theta(1)+Theta(2)*Theta(1)+Theta(2)*Phi);
 ACVF_1=Phi*ACVF_0+sigmasq*(Theta(1)+Theta(2)*Theta(1)+Theta(2));
 ACVF_2=Phi*ACVF_1+Theta(2)*sigmasq;
  ACVF=zeros(h+1,1);
  ACVF(1:3)=[ACVF_0;ACVF_1;ACVF_2];
      if h>2
          for i=4:h+1
              ACVF(i)=Phi*ACVF(i-1);
          end
      end
  ACVF=ACVF(1:h+1);
end