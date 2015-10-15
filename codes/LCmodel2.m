function EstS =  LCmodel2(yearStart,yearEnd,yearv,alpha,omega,agev,qx)


% USAGE :
% Estimates the parameters of the LC model using Poisson likelihood function

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
% EstS.model - name of model estimated
% EstS.alpha - youngest age involved (Lee and Carter are bad for fitting the
%      entire age range, so start somewhere. Say 50 or 65.
% EstS.omega - oldest age could be 85 or 90
% EstS.yearStart - first year for which estimation should be performed
% EstS.yearEnd - last year for which estimation should be made
% EstS.nbYrEst - returns the number of years over which the estimation was done
% EstS.nbAges - number of ages for which the esimtation is made
% EstS.tEst - time required for the estimation
% EstS.axhat - Beta 1 parameter
% EstS.bxhat - Beta 2 parameter
% EstS.kthat - Kappa parameter
% EstS.resids - residuals defined as difference between actual mortality and
% estimated mortality for the block of years and ages considered

retidx=find(agev==alpha);
oldidx=find(agev==omega);

ySidx=find(yearv==yearStart); %yearly Starting index
yEidx=find(yearv==yearEnd); % yearly Ending index

qxM = qx(retidx:oldidx,ySidx:yEidx);

% If there is a zero in the matrix can lead sometimes to problems
e=qxM<0.0000001; qxM(e)=0.00001;

[ExtTrue,DxtTrue] = getNbDeath(qxM); % construct lx and actual death for all ages

[nbAges,nbYrEst] = size(qxM);


% define a structure with general parameters
EstS.model = 'log(qxt)=a_x + b_x k_t Poisson likelihood estimation';
EstS.alpha = alpha;
EstS.omega = omega;
EstS.yearStart = yearStart;
EstS.yearEnd = yearEnd;
EstS.nbYrEst = nbYrEst;
EstS.nbAges  = nbAges;


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' get the Lee-Carter 1992 SVD decomposition as starting value     ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

tic
[axhat,bxhat,kthat,qxMhat] = LC1(qxM);
test=toc;
fprintf('SVD took %8.4f seconds\n',test)

qxMhat = getLC1_qxMhat(axhat,bxhat,kthat,nbAges,nbYrEst);

%   disp('ahat bxhat kthat from SVD')
%   niceprint([axhat(1:10),bxhat(1:10),kthat(1:10)'])

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' implement now the Poisson estimation using Goodman etc    ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

xSVD = [axhat(:); bxhat(:); kthat(:)];

NbParam=length(xSVD);

tic
fprintf('Compute Poisson LL using Goodman\n')
x = xSVD;
dist=10^10;iterCtr=0;
while dist>10^(-10) && iterCtr<10000
    iterCtr=iterCtr+1;
    [x1,axhat1,bxhat1,kthat1] = LCGoodmanPoiss(x,DxtTrue,ExtTrue,nbAges,nbYrEst);
    % as first described in Brouhns, Denuit, Vermunt
    dist = sum(abs(x1-x));
    x=x1;
end
EstS.tEst=toc;
fprintf('PoissonLL took %8.4f seconds\n',EstS.tEst)
EstS.dist=dist;
EstS.iterCtr=iterCtr;
EstS.axhat=axhat1;
EstS.bxhat=bxhat1;
EstS.kthat=kthat1;

qxMhat1 = getLC1_qxMhat(axhat1,bxhat1,kthat1,nbAges,nbYrEst);

%   disp('ahat bxhat kthat from LC Poisson iteration')
%   niceprint([axhat1(1:10),bxhat1(1:10),kthat1(1:10)'])

resids = qxM-qxMhat1;
EstS.resids=resids;
