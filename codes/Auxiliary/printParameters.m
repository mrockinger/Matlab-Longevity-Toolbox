function  printParameters(axhat,bxhat,kthat,age1,yearv,qxM)

nbYrEst = length(yearv);


h=figure();
subplot(2,2,1);
for t=1:nbYrEst
  plot(age1(:,1),qxM(:,t),'+k'); hold on
end
plot(age1(:,1),exp(axhat),'+r'); title('exp(axhat)'); hold off
title('\beta_x^1 hat')


subplot(2,2,2);
plot(bxhat);
title('\beta_x^2 hat')


subplot(2,2,3);
plot(yearv,kthat,'or'); hold on
plot(yearv,kthat,'xk'); hold off
title('\kappa_t hat')

print(h,'-dpdf','LC1Parameters.pdf');
