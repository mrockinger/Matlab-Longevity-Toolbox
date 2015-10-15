function x=ilogit(y)
% computes the inverse of logit
ey=exp(y);
x=ey./(1+ey);