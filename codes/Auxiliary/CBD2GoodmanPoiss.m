function x1=CBD2GoodmanPoiss(x,Dxt,Ext,N,T,age,Nc,Ns1,Ns2)
    
% Provides the estimated updated parameters via Goodman type iterations
% USAGE:
%   [kappa1_new(:); kappa2_new(:); gz_new(:)] =CBD2GoodmanPoiss(x,Dxt,Ext,N,T,age,Nc,Ns1,Ns2)
%
% INPUTS:
%
% x,Dxt,Ext,N,T,age,Nc,Ns1,Ns2

% OUTPUTS:
%
%   kappa1_new(:) = Kappa_1 updated parameter
%   kappa2_new(:) = Kappa_2 updated parameter
%   gz_new(:) = Gamma_3 updated parameter
%
% COMMENTS:
%   The CBD2 model estimated is logit q_x,t = kapa(1)_t + kappa(2)_t(x-xbar) +
%   gamma(3)_t-x + eps_x,t
%
%
k1 = x(1:T)';
k2 = x(T+1:2*T)';
gz = x(2*T+1:2*T+Nc);

% k1
% k2
% gz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update k1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num,den]=getDerivk1_CBD2(k1,k2,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);

k1_new = k1 - num./den;

%[Dhat,dy]=updatedydy_CBD2(k1_new,k2,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update k2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num,den]=getDerivk2_CBD2(k1_new,k2,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);

k2_new = k2 - num./den;

%[Dhat,dy]=updatedydy_CBD2(k1_new,k2_new,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update gz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num,den]=getDerivgz_CBD2(k1_new,k2_new,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);

gz_new = gz - num./den;

x1 = [ k1_new(:); k2_new(:); gz_new(:) ];

end