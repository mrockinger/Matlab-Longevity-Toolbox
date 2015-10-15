function qxMhat = getLC1_qxMhat(axhat,bxhat,kthat,N,T )
% computes mortality for parameters coming out of Lee-Carter estimation

eh = ones(1,T);
ev = ones(N,1);

qxMhat = exp( axhat(:,eh)+bxhat(:,eh).*kthat(ev,:) );