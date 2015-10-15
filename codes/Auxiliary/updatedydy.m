function [yh1,dy]=updatedydy(Dxt,Ext,N,Nc,Ns1,Ns2,T,ax,bx,kt,bx2,gz)
% computes number of death from given parameters, ax,bx,k_t,bx2, gam_z
% as well as differential of deaths 
eh = ones(1,T);
ev = ones(N,1);

X1    =  ax(:,eh) + bx(:,eh).*kt(ev,:) ;

X2 = zeros(N,T);

for cohidx=1:Nc  
  
  if cohidx<=N-Ns1
    z=N-Ns1-cohidx;
    for ctr=1:min(T,N-z)
      %X2(z+ctr,ctr) = cohidx %for control purposes
      X2(z+ctr,ctr) = bx2(z+ctr)*gz(cohidx); % x is a scalar !!
    end
  else
    z=cohidx-(N-Ns2)+1;
    for ctr=1:min(N,T-z+1)
      %X2(ctr,z+ctr-1) = cohidx %for control purposes
      X2(ctr,z+ctr-1) = bx2(ctr)*gz(cohidx);
    end
  end
  
end

yh1 = exp( X1 + X2 ).*Ext;
dy  = Dxt - yh1; % corresponds to y - y_hat

