function ACVF=ACVF_AR1(h,Phi,sigmasq)
% This computes the ACVF at lag h for a AR(1) model and returns a
% vector of all lags up to h.

if ~(length(Phi)==1)
    error('Length of coefficient vector is not 1')
else
   
  ACVF_0=(sigmasq)/(1-Phi^2);
  ACVF=zeros(h+1,1);
  ACVF(1)=[ACVF_0];
     if h>0
          for i=2:h+1
              ACVF(i)=Phi^(i-1)*ACVF_0;
          end
      end
  ACVF=ACVF(1:h+1);
end