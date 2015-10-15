function [kappa1v,kappa2v,kappa3v]=getInitKappasq(Ext,Dxt,age,nbAges,nbYrEst)
% computes starting values for the model with three kappas

m   = Dxt./Ext; 
q   = 1-exp(-m);
LHS = log( q./(1-q) );

kappa1v = zeros(nbYrEst,1);
kappa2v = zeros(nbYrEst,1);
kappa3v = zeros(nbYrEst,1);

agec  = age-mean(age);
agesq = agec.^2;
agesq = agesq - mean(agesq);

for t=1:nbYrEst
  
  y = LHS(:,t);
  x = [ones(nbAges,1),agec,agesq];
  res=ols(y,x);
  
  kappa1v(t)=res.beta(1);
  kappa2v(t)=res.beta(2);
  kappa3v(t)=res.beta(3);
  
end
