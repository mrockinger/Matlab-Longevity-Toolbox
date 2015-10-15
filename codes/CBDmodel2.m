function EstS =  CBDmodel2(yearStart,yearEnd,yearv,alpha,omega,agev,qx)

% logit(mortality)=kappa1_t + kappa2_tx(age-mean )+gamma_z

% USAGE :
% estimate Cairns, Blake and Dowd model with cohort effects using OLS, then with iterative scheme of Goodman

% INPUTS :
% yearStart - Year of start of the estimation of the parameters
% yearEnd - Year of end of the estimation of the parameters
% yearv    - vector of all possible years. Must contain yearStart and yearEnd
% alpha - Minimum age involved in estimation
% omega - Maximum age allowed by the model
% agev - vector of all possible ages
% qx - a matrix with as many rows as there are ages (in agev) and as many columns as there
% are years (in yearv)
%
%OUTPUTS :
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
% EstS.gz3 - Gamma parameter corresponds to a cohort effect
% EstS.resids - residuals defined as difference between actual mortality and
%       estimated mortality for the block of years and ages considered


Ns1 = 5; % skip Ns1 cohorts at the beginning of sample
Ns2 = 5; % skip Ns2 cohorts from the end of sample


retidx=find(agev==alpha);
oldidx=find(agev==omega);

ySidx=find(yearv==yearStart); %yearly Starting index
yEidx=find(yearv==yearEnd); % yearly Ending index

agec   = agev(retidx:oldidx); % conform age and year with dimension of q matrix
yearc  = yearv(ySidx:yEidx);

qxM = qx(retidx:oldidx,ySidx:yEidx);
[nbAges,nbYrEst] = size(qxM);

% If there is a zero in the matrix can lead sometimes to problems
e=qxM<0.0000001; qxM(e)=0.00001;

[ExtTrue,DxtTrue] = getNbDeath(qxM); % from given qxM construct lx and actual death

% define a structure with general parameters
EstS.model = 'logit(q_{xt})=\kappa^1(t)+\kappa^2(t)(x-bar x)+gam_z';
EstS.alpha = alpha;
EstS.omega = omega;
EstS.yearStart = yearStart;
EstS.yearEnd = yearEnd;
EstS.Ns1 = Ns1;
EstS.Ns2 = Ns2;
EstS.nbYr = nbYrEst;
EstS.nbAges  = nbAges;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' estimate CBD model: logt(q_{xt})=\kappa^1(t)+\kappa^2(t)(x-\bar x)    ')
disp(' get OLS starting values    ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

[k1,k2] = getInitKappa(ExtTrue,DxtTrue,agec,nbAges,nbYrEst); % here perform OLS estimation

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' estimate CBD model: logt(q_{xt})=\kappa^1(t)+\kappa^2(t)(x-\bar x)    ')
disp(' perform iterations à la Goodman of log-likelihood.    ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


x0=[k1(:);k2(:)];

dist=10^10;iterCtr=0; maxiter=10000;
while dist>0.0001 && iterCtr<maxiter
    
    iterCtr=iterCtr+1;
    
    x1 = CBD1GoodmanPoiss(x0,DxtTrue,ExtTrue,nbAges,nbYrEst,agec);
    
    dist = sum(abs(x1-x0));
    x0=x1;
end
%fprintf('sum of absolute deviations of mortality %8.4f\n',dist)
%fprintf('counter of how many iteration needed    %8.4f\n',iterCtr)
if iterCtr==maxiter; warning('In CBDmodel2 in first model estimation hit maxiter'); end

k1Good = x1(1:nbYrEst)';
k2Good = x1(nbYrEst+1:2*nbYrEst)';

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' estimate CBD model: logt(q_{xt})=\kappa^1(t)+\kappa^2(t)(x-bar x)+gam_z   ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
tic

% number of cohorts:
nbcohorts = nbAges-Ns1+nbYrEst-Ns2-1;
nbcohorts

x0=[k1Good(:);k2Good(:);0.1*ones(nbcohorts,1)];

tic;
dist=10^10;iterCtr=0;
while dist>0.0000001 && iterCtr<40000
    
    iterCtr=iterCtr+1;
    
    x1 = CBD2GoodmanPoiss(x0,DxtTrue,ExtTrue,nbAges,nbYrEst,agec,nbcohorts,Ns1,Ns2);
    
    dist = sum(abs(x1-x0));
    x0=x1;
end
EstS.tEst=toc;
EstS.dist=dist;
EstS.iterCtr=iterCtr;

k1  = x1(1:nbYrEst)';
k2  = x1(nbYrEst+1:2*nbYrEst)';
gz3 = x1(2*nbYrEst+1:2*nbYrEst+nbcohorts);

% now standardize the cohort effect
[k1,k2,gz3] = stdize_CBD2(k1,k2,gz3,agec,yearv);



EstS.k1=k1;
EstS.k2=k2;
EstS.gz3=gz3;

[~,~,qxMhat,~]=updatedydy_CBD2(k1,k2,gz3,DxtTrue,ExtTrue,...
    nbAges,nbYrEst,agec,nbcohorts,Ns1,Ns2);

resids = qxM - qxMhat;
EstS.resids=resids;




