function x1=CBD3GoodmanPoiss(x,Dxt,Ext,N,T,age,Nc,Ns1,Ns2)
% estimates the updated parameters via Goodman type iterations

k1 = x(1:T)';
k2 = x(T+1:2*T)';
k3 = x(2*T+1:3*T)';
gz = x(3*T+1:3*T+Nc);

% k1
% k2
% gz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update k1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num,den]=getDerivk1_CBD3(k1,k2,k3,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);

k1_new = k1 - num./den;

%[Dhat,dy]=updatedydy_CBD3(k1_new,k2,k3,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update k2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num,den]=getDerivk2_CBD3(k1_new,k2,k3,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);

k2_new = k2 - num./den;

%[Dhat,dy]=updatedydy_CBD3(k1_new,k2_new,k3,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update k3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num,den]=getDerivk3_CBD3(k1_new,k2_new,k3,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);

k3_new = k3 - num./den;

%[Dhat,dy]=updatedydy_CBD3(k1_new,k2_new,k3,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update gz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num,den]=getDerivgz_CBD3(k1_new,k2_new,k3_new,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);

gz_new = gz - num./den;





x1 = [ k1_new(:); k2_new(:); k3_new(:); gz_new(:) ];