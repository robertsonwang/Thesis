function ACVF=ACVF_ARMA11(h,Phi,Theta,sigmasq)
% This computes the ACVF at lag h for a ARMA(1,1) model and returns a
% vector of all lags up to h.

if ~(length(Phi)==1) 
    error('Length of coefficient vector is not 1')
end
if ~(length(Theta)==1)
    error('Length of coefficient vector is not 1')
end


  ACVF_0=(sigmasq)*(1+((Theta+Phi)^2)/(1-Phi^2));
  ACVF_1=(sigmasq)*(Theta+Phi+(((Theta+Phi)^2)*Phi)/(1-Phi^2));
  ACVF=zeros(h+1,1);
  ACVF(1:2)=[ACVF_0;ACVF_1];
      if h>1
          for i=3:h+1
              ACVF(i)=(Phi^(i-2))*ACVF_1;
          end
      end
  ACVF=ACVF(1:h+1);
end