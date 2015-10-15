function [x,ax_new,bx_new,kt_new]=LCGoodmanPoiss(x,Dxt,Ext,Nab,Nk)
% assumes that x contains Nab elements ax = Nab elements bx and Nk kappas
% estimates the ax, bx and kt elements via iteration
% likelihood set for Poisson density. 
ax = x(1:Nab);
bx = x(Nab+1:2*Nab);
kt = x(2*Nab+1:2*Nab+Nk)';

%disp('in function. ax, bx, kt')
%niceprint([ax(1:10),bx(1:10),kt(1:10)'])
%stop


eh=ones(1,Nk);
ev=ones(Nab,1);

mxt = exp( ax(:,eh)+bx(:,eh).*kt(ev,:) );

%niceprint(Ext(1:5,1:5))

Dxhat0 = Ext.*mxt;
% disp('dxhat')
% niceprint(Dxhat0(1:5,1:5))
% stop

ax_new = ax + sum(Dxt-Dxhat0,2)./sum(Dxhat0,2);

mxt = exp( ax_new(:,eh)+bx(:,eh).*kt(ev,:) );

Dxhat1 = Ext.*mxt;

kt_new = kt + sum((Dxt-Dxhat1).*bx(:,eh))./sum(Dxhat0.*(bx(:,eh).^2));

kt_new = kt_new - mean(kt_new(:));

mxt = exp( ax_new(:,eh)+bx(:,eh).*kt_new(ev,:) );

Dxhat2 = Ext.*mxt;

bx_new = bx + sum((Dxt-Dxhat2).*kt_new(ev,:),2)./sum(Dxhat2.*(kt_new(ev,:)).^2,2); 


% standardization procedure
sumb  = sum(bx_new);
bx_new = bx_new/sumb; kt_new = kt_new*sumb;
sumk  = mean(kt_new,2); % this is a negligible number 10^-14
kt_new = kt_new-sumk;
ax_new = ax_new + bx_new*sumk;

x=[ax_new(:); bx_new(:); kt_new(:)];
