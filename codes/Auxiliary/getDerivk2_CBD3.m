function [num,den]=getDerivk2_CBD3(k1,k2,k3,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2)

[Dhat,dy,q,m,f]=updatedydy_CBD3(k1,k2,k3,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);

x2=age-mean(age);
x2=x2(:,ones(1,T));

X = Ext.*dy.*q.*x2./Dhat;
num=sum(X);
% 
% h=0.0000001;
% 
% [Dhatp,dyp,qp,mp]=updatedydy_CBD3(k1,k2+h,k3,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);
% Xp=Ext.*dyp.*qp.*x2./Dhatp;
% 
% [Dhatn,dyn,qn,mn]=updatedydy_CBD3(k1,k2-h,k3,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);
% Xn=Ext.*dyn.*qn.*x2./Dhatn;
% 
% den=(sum(Xp)-sum(Xn))/(2*h);
% 
% disp('Den computed numerically')
% den


%tic
% now compute den via exact second order derivative
X1 = -Dxt.*Ext.*(q.^2)./(Dhat.^2);
X2 = dy.*q./(Dhat.*(1+exp(f)));
% X3 is 0 since second derivative of 1 is 0
den=sum(Ext.*((X1+X2).*(x2).^2));

%disp('den computed analytically')
%den
