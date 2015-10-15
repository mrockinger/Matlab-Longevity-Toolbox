function [Dhat,dy,q,m]=updatedydy_CBD1(k1,k2,Dxt,Ext,N,T,age)
% computes number of death from given parameters, k1, k2
% as well as differential of deaths 

k1=k1(:)';  k2=k2(:)'; % make sure these are row vectors

eh = ones(1,T);
ev = ones(N,1);

agec = age-mean(age);

f = k1(ev,:) + k2(ev,:).*agec(:,eh);

q = exp(f)./(1+exp(f));

m = -log( 1-q );

Dhat = Ext.*m;

dy  = Dxt - Dhat; % corresponds to y - y_hat

