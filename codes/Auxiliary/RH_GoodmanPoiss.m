function x=RH_GoodmanPoiss(x,Dxt,Ext,N,Nc,Ns1,Ns2,T)
% assumes that x contains
% Nb number of elements in bx
% Nc number of cohorts. The first cohort started in the past. One needs to
% know where this was in the matrix. Ns is the number of lines one must
% skip from end of mortalities to get the first cohort in end-Ns. There are
% as many iotas as there are Nc.
% Nk number of kappas
% Dxt is the number of actual deceased

% estimates the b0 (coeff of ioata), b1, and kt elements via iteration.

% [axhat1(:);bxhat1(:);kthat1(:);bxhat2(:);gzhat(:)];

p1=1;    p2=N;     ax  = x(p1:p2);
p1=p2+1; p2=2*N;   bx  = x(p1:p2);
p1=p2+1; p2=p2+T;  kt  = x(p1:p2);
p1=p2+1; p2=p2+N;  bx2 = x(p1:p2);
p1=p2+1; p2=p2+Nc; gz  = x(p1:p2);

kt=kt'; % this should be a row vector

% 
% figure()
% subplot(3,2,1)
% plot(ax)
% title('ax')
% subplot(3,2,2)
% plot(bx)
% title('bx')
% subplot(3,2,3)
% plot(kt)
% title('kt')
% subplot(3,2,4)
% plot(bx2)
% title('bx2')
% subplot(3,2,5)
% plot(gz)
% title('gamma_z')


eh = ones(1,T);
ev = ones(N,1);

% get estimate of number of death
[Dxhat,dy] = updatedydy(Dxt,Ext,N,Nc,Ns1,Ns2,T,ax,bx,kt,bx2,gz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update ax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax_new = ax + sum(dy,2)./sum(Dxhat,2);

[Dxhat,dy] = updatedydy(Dxt,Ext,N,Nc,Ns1,Ns2,T,ax_new,bx,kt,bx2,gz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update kt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kt_new = kt + sum(dy.*bx(:,eh))./sum(Dxhat.*(bx(:,eh).^2));

[Dxhat,dy] = updatedydy(Dxt,Ext,N,Nc,Ns1,Ns2,T,ax_new,bx,kt_new,bx2,gz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update bx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bx_new = bx + sum(dy.*kt_new(ev,:),2)./sum(Dxhat.*(kt_new(ev,:)).^2,2);

% rescale bx and kappa
kt_new = kt_new-mean(kt_new,2); % need to recenter immediately
c = sum(bx_new);
bx_new = bx_new/c;
kt_new = kt_new*c;

[Dxhat,dy] = updatedydy(Dxt,Ext,N,Nc,Ns1,Ns2,T,ax_new,bx_new,kt_new,bx2,gz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update gamaz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dxhat = Dxhat/10000;
% dy    = dy/10000;

Nupdate = dy  .* bx2(:,eh); % stuff to update numerator
Dupdate = Dxhat .* (bx2(:,eh).^2); % stuff to update denominator

SNupdate=zeros(Nc,1);
SDupdate=zeros(Nc,1);

for cohidx=1:Nc
  
  if cohidx<=N-Ns1
    z=N-Ns1-cohidx;
    for ctr=1:min(T,N-z)
      %X2(z+ctr,ctr) = cohidx %for control purposes
      SNupdate(cohidx) = SNupdate(cohidx) + Nupdate(z+ctr,ctr);
      SDupdate(cohidx) = SDupdate(cohidx) + Dupdate(z+ctr,ctr);
    end
  else
    z=cohidx-(N-Ns2)+1;
    for ctr=1:min(N,T-z+1)
      %X2(ctr,z+ctr-1) = cohidx %for control purposes
      SNupdate(cohidx) = SNupdate(cohidx) + Nupdate(ctr,z+ctr-1);
      SDupdate(cohidx) = SDupdate(cohidx) + Dupdate(ctr,z+ctr-1);
    end
  end
  
end

%niceprint([SNupdate,SDupdate])

gz_new = gz + SNupdate./SDupdate;
gz_new = gz_new-mean(gz_new);
% figure()
% plot(gz_new)
% title('gz')

[Dxhat,dy] = updatedydy(Dxt,Ext,N,Nc,Ns1,Ns2,T,ax_new,bx_new,kt_new,bx2,gz_new);

% Dxhat = Dxhat/10000;
% dy    = dy/10000;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update bx2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X2N = zeros(N,T);
X2D = zeros(N,T);

for cohidx=1:Nc  
  
  if cohidx<=N-Ns1
    z=N-Ns1-cohidx;
    for ctr=1:min(T,N-z)
      %X2(z+ctr,ctr) = cohidx %for control purposes
      X2N(z+ctr,ctr) = gz_new(cohidx); % x is a scalar !!
      X2D(z+ctr,ctr) = gz_new(cohidx)^2; % x is a scalar !!      
    end
  else
    z=cohidx-(N-Ns2)+1;
    for ctr=1:min(N,T-z+1)
      %X2(ctr,z+ctr-1) = cohidx %for control purposes
      X2N(ctr,z+ctr-1) = gz_new(cohidx);
      X2D(ctr,z+ctr-1) = gz_new(cohidx)^2; % x is a scalar !!        
    end
  end
  
end

num = sum(dy.*X2N,2);
dum = sum(Dxhat.*X2D,2);
%niceprint([num,dum])
bx2_new = bx2 + num./dum;

% figure()
% plot(bx2_new)
% title('bx2_new')

 gz_new = gz_new - mean(gz_new); % need to recenter immediately
 c = sum(bx2_new);
 bx2_new = bx2_new/c;
 gz_new = gz_new*c;


x = [ax_new(:);bx_new(:);kt_new(:);bx2_new(:);gz_new(:)];

% 
% figure()
% subplot(3,2,1)
% plot(ax_new)
% title('ax')
% subplot(3,2,2)
% plot(bx_new)
% title('bx')
% subplot(3,2,3)
% plot(kt_new)
% title('kt')
% subplot(3,2,4)
% plot(bx2_new)
% title('bx2')
% subplot(3,2,5)
% plot(gz_new)
% title('gamma_z')


