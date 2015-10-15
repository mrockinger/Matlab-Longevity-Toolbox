function plotRHestimates(x,N,Nc,Ns1,Ns2,T)
% extracts from x0 and x1 the ax, bx etc components and plot the variosu
% estimates


p1=1;    p2=N;     ax  = x(p1:p2);
p1=p2+1; p2=2*N;   bx  = x(p1:p2);
p1=p2+1; p2=p2+T;  kt  = x(p1:p2);
p1=p2+1; p2=p2+N;  bx2 = x(p1:p2);
p1=p2+1; p2=p2+Nc; gz  = x(p1:p2);

h=figure();
subplot(3,2,1)
plot(ax)
title('ax')
subplot(3,2,2)
plot(bx)
title('bx')
subplot(3,2,3)
plot(kt)
title('kt')
subplot(3,2,4)
plot(bx2)
title('bx2')
subplot(3,2,5)
plot(gz)
title('gz')

