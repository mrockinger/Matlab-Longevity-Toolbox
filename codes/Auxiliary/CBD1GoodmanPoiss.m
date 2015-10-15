function x1=CBD1GoodmanPoiss(x,Dxt,Ext,N,T,age)

% Provides the estimated updated parameters via Goodman type iterations
% USAGE:
%   [kappa1_new(:); kappa2_new(:)] =CBD1GoodmanPoiss(x,Dxt,Ext,N,T,age)
%
% INPUTS:
%
% x,Dxt,Ext,N,T,age
%
% OUTPUTS:
%
%   kappa1_new(:) = Kappa_1 updated parameter
%   kappa2_new(:) = Kappa_2 updated parameter
%
% COMMENTS:
%   The CBD1 model estimated is logit q_x,t = kapa(1)_t + kappa(2)_t(x-xbar) +
%   eps_x,t
%

k1 = x(1:T)';
k2 = x(T+1:2*T)';

%[Dhat,dy] = updatedydy_CBD1(k1,k2,Dxt,Ext,N,T,age);

% eh = ones(1,T);
% ev = ones(N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update k1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num,den]=getDerivk1_CBD1(k1,k2,Dxt,Ext,N,T,age);

k1_new = k1 - num./den;

[Dhat,dy]=updatedydy_CBD1(k1_new,k2,Dxt,Ext,N,T,age);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update k2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num,den]=getDerivk2_CBD1(k1_new,k2,Dxt,Ext,N,T,age);

k2_new = k2 - num./den;

x1 = [ k1_new(:); k2_new(:) ];