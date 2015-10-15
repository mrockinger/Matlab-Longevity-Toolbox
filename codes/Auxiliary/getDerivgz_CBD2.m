function [num,den]=getDerivgz_CBD2(k1,k2,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2)

% compute partial derivative of Deviance wrt gz
[Dhat,dy,q,m,f]=updatedydy_CBD2(k1,k2,gz,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);

X=Ext.*dy.*q./Dhat;

num = zeros(Nc,1); % stuff to update numerator

for cohidx=1:Nc
  
  if cohidx<=N-Ns1
    z=N-Ns1-cohidx;
    for ctr=1:min(T,N-z)
      %X2(z+ctr,ctr) = cohidx %for control purposes
      num(cohidx) = num(cohidx) + X(z+ctr,ctr);

    end
  else
    z=cohidx-(N-Ns2)+1;
    for ctr=1:min(N,T-z+1)
      %X2(ctr,z+ctr-1) = cohidx %for control purposes
      num(cohidx) = num(cohidx) + X(ctr,z+ctr-1);      
    end
  end
  
end

%tic
% now compute den via exact second order derivative
% compute the second partial derivative for denuminator
X1 = -Dxt.*Ext.*(q.^2)./(Dhat.^2);
X2 = dy.*q./(Dhat.*(1+exp(f)));
%X3 = dy.*q./Dhat;
% X3 is 0 since second derivative of 1 is 0
X=(Ext.*(X1+X2));

den = zeros(Nc,1); % stuff to update numerator

for cohidx=1:Nc
  
  if cohidx<=N-Ns1
    z=N-Ns1-cohidx;
    for ctr=1:min(T,N-z)
      %X2(z+ctr,ctr) = cohidx %for control purposes
      den(cohidx) = den(cohidx) + X(z+ctr,ctr);

    end
  else
    z=cohidx-(N-Ns2)+1;
    for ctr=1:min(N,T-z+1)
      %X2(ctr,z+ctr-1) = cohidx %for control purposes
      den(cohidx) = den(cohidx) + X(ctr,z+ctr-1);      
    end
  end
  
end


% 
% 
% 
% disp('den computed analytically')
% den'
% toc
% 
% 
% tic
% 
% 
% % compute numerically second partial derivative of Deviance wrt gz
% h=0.0000001;
% 
% [Dhatp,dyp,qp,mp]=updatedydy_CBD2(k1,k2,gz+h,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);
% Xp=Ext.*dyp.*qp./Dhatp;
% [Dhatn,dyn,qn,mn]=updatedydy_CBD2(k1,k2,gz-h,Dxt,Ext,N,T,age,Nc,Ns1,Ns2);
% Xn=Ext.*dyn.*qn./Dhatn;
% 
% denp = zeros(Nc,1); % stuff to update numerator
% denn = zeros(Nc,1); % stuff to update numerator
% 
% for cohidx=1:Nc
%   
%   if cohidx<=N-Ns1
%     z=N-Ns1-cohidx;
%     for ctr=1:min(T,N-z)
%       %X2(z+ctr,ctr) = cohidx %for control purposes
%       denp(cohidx) = denp(cohidx) + Xp(z+ctr,ctr);
%       denn(cohidx) = denn(cohidx) + Xn(z+ctr,ctr);
%     end
%   else
%     z=cohidx-(N-Ns2)+1;
%     for ctr=1:min(N,T-z+1)
%       %X2(ctr,z+ctr-1) = cohidx %for control purposes
%       denp(cohidx) = denp(cohidx) + Xp(ctr,z+ctr-1);      
%       denn(cohidx) = denn(cohidx) + Xn(ctr,z+ctr-1);            
%     end
%   end
%   
% end
% 
% 
% den=(denp-denn)/(2*h);
% 
% disp('den computed numerically')
% den'
% toc
% stop
