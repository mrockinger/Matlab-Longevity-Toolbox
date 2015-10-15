function [Ext,Dxt] = getNbDeath(qxM)
% estimate the number of death for each age/year as well as the number of
% exposures (that is the population at beginning of year). Start with
% cohorts of 100'000 in the first year
[nbAges,nbYrEst]=size(qxM);

Ext=zeros(nbAges+1,nbYrEst);

Ext(1,:) = 100000;
for t=2:nbAges+1
  Ext(t,:) = Ext(t-1,:).*(1-qxM(t-1,:));
end
Dxt=Ext(1:end-1,:)-Ext(2:end,:); % number of death
Ext=Ext(1:end-1,:); % exposed to risk