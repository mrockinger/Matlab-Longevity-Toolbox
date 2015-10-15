function [num,den]=getDerivk1_CBD1(k1,k2,Dxt,Ext,N,T,age)

[Dhat,dy,q,m]=updatedydy_CBD1(k1,k2,Dxt,Ext,N,T,age);

X=Ext.*dy.*q./Dhat;
num=sum(X);

h=0.0000001;

[Dhatp,dyp,qp,mp]=updatedydy_CBD1(k1+h,k2,Dxt,Ext,N,T,age);
Xp=Ext.*dyp.*qp./Dhatp;
[Dhatn,dyn,qn,mn]=updatedydy_CBD1(k1-h,k2,Dxt,Ext,N,T,age);
Xn=Ext.*dyn.*qn./Dhatn;
den=(sum(Xp)-sum(Xn))/(2*h);




