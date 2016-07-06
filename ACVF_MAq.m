function ACVF=ACVF_MAq(q,h,Theta,sigmasq)
% This computes the ACVF at lag h for a MA(q) model and returns a
% vector of all lags up to h.

if ~(length(Theta)==q)
    error('Length of coefficient vector is not correct')
else
   
    
       s=1;     
    for j=1:q
        s=s+Theta(j)^2;
    end
    ACVF_0=sigmasq*s;
    ACVF=zeros(h+1,1);
    ACVF(1)=ACVF_0;
          for i=2:h+1
              if (i-1)<=q
                  sum=Theta(i-1);
                  for k=i:q
                      sum=sum+Theta(k)*Theta(k-i+1);
                  end
              ACVF(i)=sigmasq*sum;
              else
              ACVF(i)=0;
              end
          end
      end
  ACVF=ACVF(1:h+1);
end