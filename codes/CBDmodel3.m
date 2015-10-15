function EstS =  CBDmodel2(yearStart,yearEnd,yearv,alpha,omega,agev,qx)

% logit(mortality)=kappa1_t + kappa2_tx(age-mean ) + kappa2_tx((age-mean)^2-s) + gamma_z

% USAGE :
% estimate Cairns, Blake and Dowd model with cohort effects and quadratic term using
% OLS (to get starting values), then with iterative scheme of Goodman


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

% k1- Kappa 1 parameter
% k2 - Kappa 2 parameter
% k3 - Kappa 3 parameter
% gz3 - Gamma parameter
% resids - residuals defined as difference between actual mortality and
% estimated mortality for the block of years and ages considered

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
EstS.model = 'logit(q_{xt})=\k1(t)+\k^2(t)(x-\bar x) +\k^3(t) [(x-\bar x)^2 -sigma^2]+gz';
EstS.alpha = alpha;
EstS.omega = omega;
EstS.yearStart = yearStart;
EstS.yearEnd = yearEnd;
EstS.Ns1 = Ns1;
EstS.Ns2 = Ns2;
EstS.nbYr = nbYrEst;
EstS.nbAges  = nbAges;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' estimate CBD model: logt(q_{xt})=\kappa^1(t)+\kappa^2(t)(x-\bar x)')
disp('     +\kappa^3(t) [(x-\bar x)^2 -sigma^2] ')
disp(' get OLS starting values    ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

[k1,k2,k3]=getInitKappasq(ExtTrue,DxtTrue,agec,nbAges,nbYrEst); % here perform OLS estimation
% niceprint([[k1,k2,k3]])

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' estimate CBD model: logt(q_{xt})=\kappa^1(t)+\kappa^2(t)(x-bar x)   ')
disp('                          +\kappa^3(t) [(x-\bar x)^2 -sigma^2] +gam_z ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
tic


Ns1 = 10; % skip Ns1 cohorts at the beginning of sample
Ns2 = 10; % skip Ns2 cohorts from the end of sample

% number of cohorts:
nbcohorts = nbAges-Ns1+nbYrEst-Ns2-1;
nbcohorts

%gz_0= startvalgz_CBD2(k1Good,k2Good,DxtTrue,ExtTrue,nbAges,nbYr,agev,nbcohorts,Ns1,Ns2);

%load gzCBD2 gzCBD2;

x0=[k1(:);k2(:);k3(:);0.1*ones(nbcohorts,1)];

NbParam=length(x0);

tic;
dist=10^10; iterCtr=0;
while dist>0.00001 && iterCtr<40000
    if mod(iterCtr,5000)==0;
        fprintf('Iteration %6.0f\n',iterCtr);
    end
    iterCtr=iterCtr+1;
    
    x1 = CBD3GoodmanPoiss(x0,DxtTrue,ExtTrue,nbAges,nbYrEst,agec,nbcohorts,Ns1,Ns2);
    
    dist = sum(abs(x1-x0));
    x0=x1;
end
EstS.tEst=toc;
EstS.dist=dist;
EstS.iterCtr=iterCtr;

k1 = x1(1:nbYrEst)';
k2 = x1(nbYrEst+1:2*nbYrEst)';
k3 = x1(2*nbYrEst+1:3*nbYrEst)';
gz3 = x1(2*nbYrEst+1:2*nbYrEst+nbcohorts);


s=(1:nbYrEst)';
figure(1)
plot(s,k1,'-k'); hold on
plot(s,k1','xb'); hold off
title('kappa 1 OLS, with Goodman type iterations and cohort model' )

figure(2)
plot(s,k2,'-k'); hold on
plot(s,k2','xb'); hold off
title('kappa 2 OLS, with Goodman type iterations and cohort model' )

figure(3)
plot(s,k3,'-k'); hold on
plot(s,k3','xb'); hold off
title('kappa 2 OLS, with Goodman type iterations and cohort model' )

s=(1:nbcohorts);
figure(4)
%plot(s,gz_0,'-k'); hold off
plot(s,gz3,'xb'); hold off
title('gz cohort effect (intial value)and Poisson-MLand with Goodman type iterations')


[~,~,qxMhat,~,~]=updatedydy_CBD3(k1,k2,k3,gz3,DxtTrue,ExtTrue,...
    nbAges,nbYrEst,agec,nbcohorts,Ns1,Ns2);
plot(gz3)

[k1CBD3A,k2CBD3A,k3CBD3A,gzCBD3A] = stdize_CBD3(k1,k2,k3,gz3,agec,yearc);

[~,~,qxMhat2,~]=updatedydy_CBD3(k1CBD3A,k2CBD3A,k3CBD3A,gzCBD3A,DxtTrue,ExtTrue,...
    nbAges,nbYrEst,agec,nbcohorts,Ns1,Ns2);

niceprint(qxMhat(1:10,1:10))
disp('---')
niceprint(qxMhat2(1:10,1:10))

resids=qxM-qxMhat;
EstS.resids=resids;


