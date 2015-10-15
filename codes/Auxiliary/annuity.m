function A = annuity(qx,i)

% Annuities computed with interest rate i and given mortality rates
%
% USAGE:
%   A = annuity(qx,i)
%
% INPUTS:
%   qx  = the mortality probabilities.
%   age the corresponding ages. retage is time when annuity should be
%   computed
%   i = used mortality rates
%
% OUTPUTS:
%   A - Annuities due in a given currency.

G = cumprod(1-qx); %survival probabilities
Nbannuities = size(qx,1);

v  = 1/(1+i); % discount factors
df = v.^(1:Nbannuities); % row vector of discount factors.

A = df*G; %annuities
