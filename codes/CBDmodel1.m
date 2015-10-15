function EstS =  CBDmodel1(yearStart,yearEnd,yearv,alpha,omega,agev,qx)
% logit(q(t,x)) = \kappa_1 + \kappa2 (x-\bar x)

% USAGE :
% Estimate Cairns, Blake and Dowd model with OLS (to get starting values),
% then with iterative scheme of Goodman

% INPUTS :
% yearStart - Year of start of the estimation of the parameters
% yearEnd - Year of end of the estimation of the parameters
% yearv    - vector of all possible years. Must contain yearStart and yearEnd
% alpha - Minimum age involved in estimation
% omega - Maximum age allowed by the model
% agev - vector of all possible ages
% qx - a matrix with as many rows as there are ages (in agev) and as many columns as there
% are years (in yearv)

%OUTPUTS :
% EstS.model - name of model estimated
% EstS.alpha - youngest age involved (Lee and Carter are bad for fitting the
%      entire age range, so start somewhere. Say 50 or 65.
% EstS.omega - oldest age could be 85 or 90
% EstS.yearStart - first year for which estimation should be performed
% EstS.yearEnd - last year for which estimation should be made
% EstS.nbYrEst - returns the number of years over which the estimation was done
% EstS.nbAges - number of ages for which the esimtation is made
% EstS.Ns1 - number of cohorts to skip (bcse of small sample pb) at beginning
% EstS.Ns2 - same as Ns1 but for last cohorts
% EstS.tEst - time required for the estimation
% EstS.k1 - Kappa 1 parameter
% EstS.k2 - Kappa 2 parameter
% EstS.resids - residuals defined as difference between actual mortality and
%       estimated mortality for the block of years and ages considered

% this is only required for plots of resdiuals
Ns1 = 5; % skip Ns1 cohorts at the beginning of sample
Ns2 = 5; % skip Ns2 cohorts from the end of sample

retidx=find(agev==alpha);
oldidx=find(agev==omega);

ySidx=find(yearv==yearStart); %yearly Starting index
yEidx=find(yearv==yearEnd); % yearly Ending index

agec   = agev(retidx:oldidx); % age and year conform with dimension of q matrix
yearc  = yearv(ySidx:yEidx);

qxM = qx(retidx:oldidx,ySidx:yEidx);
[nbAges,nbYrEst] = size(qxM);

% If there is a zero in the matrix can lead sometimes to problems
e=qxM<0.0000001; qxM(e)=0.00001;

[ExtTrue,DxtTrue] = getNbDeath(qxM); % from given qxM construct lx and actual death

% define a structure with general parameters
EstS.model   = 'logit(q_{xt})=\kappa^1(t)+\kappa^2(t)(x-\bar x) ';
EstS.alpha   = alpha;
EstS.omega   = omega;
EstS.yearStart = yearStart;
EstS.yearEnd = yearEnd;
EstS.Ns1     = Ns1;
EstS.Ns2     = Ns2;
EstS.nbYrEst = nbYrEst;
EstS.nbAges  = nbAges;


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' estimate CBD model: logt(q_{xt})=\kappa^1(t)+\kappa^2(t)(x-\bar x)    ')
disp(' get OLS starting values    ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
tic

[k1_OLS,k2_OLS] = getInitKappa(ExtTrue,DxtTrue,agec,nbAges,nbYrEst); % here perform OLS estimation

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' estimate CBD model: logt(q_{xt})=\kappa^1(t)+\kappa^2(t)(x-\bar x)    ')
disp(' perform iterations à la Goodman of log-likelihood. Use numerical approx of gradient.    ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


x0=[k1_OLS(:);k2_OLS(:)];

tic;
dist=10^10;iterCtr=0;
while dist>0.000001 && iterCtr<40000
    
    iterCtr=iterCtr+1;
    
    x1 = CBD1GoodmanPoiss(x0,DxtTrue,ExtTrue,nbAges,nbYrEst,agec);
    
    dist = sum(abs(x1-x0));
    x0=x1;
end
EstS.tEst   = toc;
EstS.dist   = dist;
EstS.iterCtr= iterCtr;

k1Good = x1(1:nbYrEst)';
k2Good = x1(nbYrEst+1:2*nbYrEst)';

EstS.k1=k1Good;
EstS.k2=k2Good;
% it was verified that estimation via an optimizer yielded the same results

s=(1:nbYrEst)';
figure(1)
plot(s,k1_OLS','-k');
plot(s,k1Good','-r'); 
title('kappa 1 OLS and Poisson-ML with Goodman type iterations' )

figure(2)
plot(s,k2_OLS,'-k'); hold on
plot(s,k2Good','-r'); hold off
title('kappa 2 OLS and Poisson-ML with Goodman type iterations')

ev=ones(nbAges,1);
eh=ones(1,nbYrEst);
agecent=agec-mean(agec);
k1Good=k1Good(:)'; k2Good=k2Good(:)';

X = k1Good(ev,:)+k2Good(ev,:).*agecent(:,eh) ;

%%%!!!!!!!qxMhat=1-exp(-muxthat);
qxMhat = exp(X)./(1+exp(X));

resids=qxM-qxMhat;
EstS.resids=resids;

mxthat = -log(1-qxMhat);

%corresponding mortality
mxt = -log(1-qxM);

agec   = agev(retidx:oldidx); % age and year conform with dimension of q matrix
yearc  = yearv(ySidx:yEidx);

figure(3)
plot(agec,mxt,'+k'); hold on
plot(agec,mxthat,'xr'); hold off
legend('True values','Estimated values')
title('m_xt and m_xt hat')


[~,Dxthat] = getNbDeath(mxthat);
%size(Dxthat)

figure(4)
plot(agec,DxtTrue,'+k'); hold on
plot(agec,Dxthat,'xr'); hold off
legend('True values','Estimated values')
title('Death and estimates thereoff')

