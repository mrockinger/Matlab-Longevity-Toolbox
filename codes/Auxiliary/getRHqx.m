function qxt = getRHqx(EstS)
% returns the mortality rate for RH model

axhat1=EstS.axhat1;
bxhat1=EstS.bxhat1;
kthat1=EstS.kthat1; kthat1=kthat1';
bxhat2=EstS.bxhat2;
gzhat =EstS.gzhat;

T  = EstS.nbYrEst;
N  = EstS.nbAges;
Nc = EstS.nbCohorts;
Ns1= EstS.Ns1;
Ns2= EstS.Ns2;

eh = ones(1,T);
ev = ones(N,1);

X1    =  axhat1(:,eh) + bxhat1(:,eh).*kthat1(ev,:) ;

X2 = zeros(N,T);

for cohidx=1:Nc
  
  if cohidx<=N-Ns1
    z=N-Ns1-cohidx;
    for ctr=1:min(T,N-z)
      %X2(z+ctr,ctr) = cohidx %for control purposes
      X2(z+ctr,ctr) = bxhat2(z+ctr)*gzhat(cohidx); % x is a scalar !!
    end
  else
    z=cohidx-(N-Ns2)+1;
    for ctr=1:min(N,T-z+1)
      %X2(ctr,z+ctr-1) = cohidx %for control purposes
      X2(ctr,z+ctr-1) = bxhat2(ctr)*gzhat(cohidx);
    end
  end
  
end

qxt    = exp( X1 + X2 ); % LC and RH provide directly mortality rate
