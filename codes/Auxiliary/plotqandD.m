function plotqandD(qxM,qxMhat,agev)
% print actual and forecasted mortality as well as acutal and forcasted
% death

h=figure();

subplot(2,1,1);
plot(agev,qxM,'+k'); hold on
plot(agev,qxMhat,'xr'); hold off
title('q_x and q_x hat')


subplot(2,1,2);
[~,Dxt] = getNbDeath(qxM);
[~,Dxthat] = getNbDeath(qxMhat);

plot(agev,Dxt,'+k'); hold on
plot(agev,Dxthat,'xr'); hold off
title('D_x and D_x hat')

print(h,'-dpdf','LC1qandD.pdf');
