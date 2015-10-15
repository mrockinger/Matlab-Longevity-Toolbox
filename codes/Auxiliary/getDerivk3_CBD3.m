function [num,den]=getDerivk3_CBD3(k1,k2,k3,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2)

[Dhat,dy,q,m,f]=updatedydy_CBD3(k1,k2,k3,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);

x2 = (age-mean(age)).^2;
x3 = x2 - mean(x2);

x3 = x3(:,ones(1,T));

X = Ext.*dy.*q.*x3./Dhat;
num=sum(X);

h=0.0000001;

[Dhatp,dyp,qp,mp]=updatedydy_CBD3(k1,k2,k3+h,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);
Xp=Ext.*dyp.*qp.*x3./Dhatp;

[Dhatn,dyn,qn,mn]=updatedydy_CBD3(k1,k2,k3-h,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);
Xn=Ext.*dyn.*qn.*x3./Dhatn;

den=(sum(Xp)-sum(Xn))/(2*h);
% 
% disp('Den computed numerically')
% den


%tic
% now compute den via exact second order derivative
X1 = -Dxt.*Ext.*(q.^2)./(Dhat.^2);
X2 = dy.*q./(Dhat.*(1+exp(f)));
% X3 is 0 since second derivative of 1 is 0
den=sum(Ext.*((X1+X2).*(x3).^2));

% disp('den computed analytically')
% den
% stop