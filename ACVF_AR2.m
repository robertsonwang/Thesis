function ACVF=ACVF_AR2(h,Phi,sigmasq)
% This computes the ACVF at lag h for a AR(2) model and returns a
% vector of all lags up to h.

if ~(length(Phi)==2)
    error('Length of coefficient vector is not 2')
else
  ACVF_0=(sigmasq)/(1-((Phi(1))^2/(1-Phi(2)))-((Phi(2)*Phi(1)^2)/(1-Phi(2)))-Phi(2)^2);
  ACVF_1=(Phi(1)/(1-Phi(2)))*ACVF_0;
  ACVF_2=Phi(1)*ACVF_1+Phi(2)*ACVF_1;
  ACVF=zeros(h+1,1);
  ACVF(1:3)=[ACVF_0;ACVF_1;ACVF_2];
      if h>2
          for i=4:h+1
              ACVF(i)=Phi(1)*ACVF(i-1)+Phi(2)*ACVF(i-2);
          end
      end
  ACVF=ACVF(1:h+1);
end