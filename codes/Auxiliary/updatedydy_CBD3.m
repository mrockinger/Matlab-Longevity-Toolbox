function [Dhat,dy,q,m,f]=updatedydy_CBD3(k1,k2,k3,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2)
% computes number of death from given parameters, k1, k2
% as well as differential of deaths 
% Ns2 is normally given implicitly with the numer of elements of gz

k1=k1(:)';  k2=k2(:)'; k3=k3(:)'; % make sure these are row vectors

eh = ones(1,T);
ev = ones(N,1);

agec  = age-mean(age);
agesq = agec.^2;
agesq = agesq - mean(agesq);

X2 = zeros(N,T);

for cohidx=1:Nc  
  
  if cohidx<=N-Ns1
    z=N-Ns1-cohidx;
    for ctr=1:min(T,N-z)
      %X2(z+ctr,ctr) = cohidx %for control purposes
      X2(z+ctr,ctr) = gz(cohidx); % x is a scalar !!
    end
  else
    z=cohidx-(N-Ns1)+1;
    for ctr=1:min(N,T-z+1)
      %X2(ctr,z+ctr-1) = cohidx %for control purposes
      X2(ctr,z+ctr-1) = gz(cohidx);
    end
  end
  
end


f = k1(ev,:) + k2(ev,:).*agec(:,eh) + k3(ev,:).*agesq(:,eh) + X2;

q = exp(f)./(1+exp(f));

m = -log( 1-q );

Dhat = Ext.*m;

dy  = Dxt - Dhat; % corresponds to y - y_hat

