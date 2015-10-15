function [num,den]=getDerivk1_CBD2(k1,k2,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2)



[Dhat,dy,q,m,f]=updatedydy_CBD2(k1,k2,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);

X=Ext.*dy.*q./Dhat;
num=sum(X);

%tic
% now compute den via exact second order derivative
X1 = -Dxt.*Ext.*(q.^2)./(Dhat.^2);
X2 = dy.*q./(Dhat.*(1+exp(f)));
% X3 is 0 since second derivative of 1 is 0
den=sum(Ext.*(X1+X2));





%disp('den computed analytically')
%den
% toc
% 
% 
% tic
% h=0.0000001;
% 
% [Dhatp,dyp,qp,mp]=updatedydy_CBD2(k1+h,k2,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);
% Xp=Ext.*dyp.*qp./Dhatp;
% 
% 
% [Dhatn,dyn,qn,mn]=updatedydy_CBD2(k1-h,k2,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);
% Xn=Ext.*dyn.*qn./Dhatn;
% 
% den=(sum(Xp)-sum(Xn))/(2*h);
% disp('den computed numerically')
% den
% 
% toc
% % 
% stop
% 





