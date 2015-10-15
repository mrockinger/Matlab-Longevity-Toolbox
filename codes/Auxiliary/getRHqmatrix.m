function [qxt,Dxthat] = getRHqmatrix(x,Ext,N,Nc,Ns1,Ns2,T)
% returns the mortality rate and the actual number of death for RH model

p1=1;    p2=N;     ax  = x(p1:p2);
p1=p2+1; p2=2*N;   bx  = x(p1:p2);
p1=p2+1; p2=p2+T;  kt  = x(p1:p2)';
p1=p2+1; p2=p2+N;  bx2 = x(p1:p2);
p1=p2+1; p2=p2+Nc; gz  = x(p1:p2);

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

qxt    = exp( X1 + X2 ); % LC and RH provide directly mortality rate
Dxthat = qxt.*Ext;       % number of deaths