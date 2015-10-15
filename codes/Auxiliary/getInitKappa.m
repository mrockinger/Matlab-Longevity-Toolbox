function [kappa1v,kappa2v] = getInitKappa(Ext,Dxt,age,nbAges,nbYrEst)
% computes starting values for the two kappas.
% Model is logit(q_{x,t}) = k_t^1 + k_t^2 (x- \bar x)
% compute y_{x,t}=logit(q_{x,t}) and perform OLS regression. This should
% have the same interpretation as Gaussian count regression

m   = Dxt./Ext; 
q   = 1-exp(-m);
LHS = log( q./(1-q) );

kappa1v = zeros(nbYrEst,1);
kappa2v = zeros(nbYrEst,1);

agec = age - mean(age);

for t=1:nbYrEst
  
  y = LHS(:,t);
  x = [ones(nbAges,1),agec];
  res = ols(y,x);
  kappa1v(t) = res.beta(1);
  kappa2v(t) = res.beta(2);
  
end
