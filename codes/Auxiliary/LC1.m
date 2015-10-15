function [axhat,bxhat,kthat,qxMhat] = LC1(qxt)
% plain Lee-Carter model
% estimates log(q_{x,t} = a_x + b_x \kappa_t
% x=1,...,N, t=1,...T.
% constaints_ \sum_x bx=1, \sum_t \kappa_t=0
% Input
% Matrix NxT of mortality rates q_{x,t}, x=age, t=year
% Output
% returns axhat Nx1, bxhat Nx1, \kappa_t 1xt
% qxMhat the forecasted probabilities

% average age effect
lqxt   = log(qxt);
axhat  = mean(lqxt,2);

[N,T]  = size(qxt);

ut = bsxfun(@minus, lqxt, axhat); % center log mortality

[U,S,V] = svd(ut);
bxhat   = -U(:,1)*S(1,1);
kthat   = -V(:,1);

kthat   = kthat(:)'; % now a row vector

% standardization procedure
% sum_x bx=1, sum_t kt=0
sumb  = sum(bxhat);
bxhat = bxhat/sumb; kthat = kthat*sumb;
sumk  = mean(kthat,2); % this is a negligible number 10^-14
kthat = kthat-sumk;
axhat = axhat + bxhat*sumk;

% reconstruction of parameters

qxMhat = getLC1_qxMhat(axhat,bxhat,kthat,N,T);

