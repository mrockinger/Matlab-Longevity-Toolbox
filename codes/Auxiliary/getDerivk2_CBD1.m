function [num,den]=getDerivk2_CBD1(k1,k2,Dxt,Ext,N,T,age)

[Dhat,dy,q,m]=updatedydy_CBD1(k1,k2,Dxt,Ext,N,T,age);

x2=age-mean(age);
x2=x2(:,ones(1,T));

X = Ext.*dy.*q.*x2./Dhat;
num=sum(X);

h=0.0000001;

[Dhatp,dyp,qp,mp]=updatedydy_CBD1(k1,k2+h,Dxt,Ext,N,T,age);
Xp=Ext.*dyp.*qp.*x2./Dhatp;

[Dhatn,dyn,qn,mn]=updatedydy_CBD1(k1,k2-h,Dxt,Ext,N,T,age);
Xn=Ext.*dyn.*qn.*x2./Dhatn;

den=(sum(Xp)-sum(Xn))/(2*h);