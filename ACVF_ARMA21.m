function ACVF=ACVF_ARMA21(h,Phi,Theta,sigmasq)
% This computes the ACVF at lag h for a ARMA(2,1) model and returns a
% vector of all lags up to h.

if ~(length(Phi)==2) 
    error('Length of coefficient vector is not 2')
end
if ~(length(Theta)==1)
    error('Length of coefficient vector is not 1')
end

 ACVF_0=(1/(1-((Phi(1)^2)/(1-Phi(2)))-(Phi(2)*Phi(1)^2)/(1-Phi(2))-Phi(2)^2))*(sigmasq*(1+Theta*Phi(1)+Theta^2)+(Phi(1)*Theta*sigmasq)/(1-Phi(2))+(Phi(2)*Phi(1)*Theta*sigmasq)/(1-Phi(2)));
 ACVF_1=(Phi(1)*ACVF_0+Theta*sigmasq)/(1-Phi(2));
 ACVF=zeros(h+1,1);
 ACVF(1:2)=[ACVF_0;ACVF_1];
      if h>1
          for i=3:h+1
              ACVF(i)=Phi(1)*ACVF(i-1)+Phi(2)*ACVF(i-2);
          end
      end
  ACVF=ACVF(1:h+1);
end