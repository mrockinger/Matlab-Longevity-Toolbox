function  EstS = RHmodel(yearStart,yearEnd,yearv,alpha,omega,agev,qx)

% USAGE :
% Estimates the parameters of the Renshaw-Haberman model

% INPUTS :
% yearStart - Year of start of the estimation of the parameters
% yearEnd - Year of end of the estimation of the parameters
% yearu   - vector of all possible years. Must contain yearStart and yearEnd
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
% EstS.nbAges - number of ages for which the estimation is made
% EstS.nbcohorts - number of cohorts in sample
% EstS.Ns1 - number of cohorts to skip (bcse of small sample pb) at beginning
% EstS.Ns2 - same as Ns1 but for last cohorts
% EstS.tEst - time required for the estimation
% EstS.axhat - Beta 1 parameter
% EstS.bxhat - Beta 2 parameter
% EstS.kthat - Kappa parameter
% EstS.bxhat2 - Beta 3 parameter
% EstS.gzhat - Gamma parameter
% EstS.resids - residuals defined as difference between actual mortality and
%         estimated mortality for the block of years and ages considered


Ns1 = 5; % skip Ns1 cohorts at the beginning of sample
Ns2 = 5; % skip Ns2 cohorts from the end of sample



retidx=find(agev==alpha);
oldidx=find(agev==omega);

ySidx=find(yearv==yearStart); %yearly Starting index
yEidx=find(yearv==yearEnd); % yearly Ending index

qxM = qx(retidx:oldidx,ySidx:yEidx);
[nbAges,nbYr] = size(qxM);

% If there is a zero in the matrix can lead sometimes to problems
e=qxM<0.0000001; qxM(e)=0.00001;

[ExtTrue,DxtTrue] = getNbDeath(qxM); % from given qxM construct lx and actual death

[nbAges,nbYrEst]=size(qxM);


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('get the Lee-Carter 1992 decomposition')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

tic
[axhat,bxhat,kthat,qxMhat] = LC1(qxM);
test=toc;
fprintf('SVD took %8.4f seconds\n',test)

agec   = agev(retidx:oldidx); % age and year conform with dimension of q matrix
yearc  = yearv(ySidx:yEidx);  % 

figure(1)
for t=1:nbYr
    plot(agec,qxM(:,t),'+k'); hold on
end
plot(agec,exp(axhat),'+r'); title('exp(axhat)'); hold off
title('ax hat')

qxMhat = getLC1_qxMhat(axhat,bxhat,kthat,nbAges,nbYrEst);


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' implement now the Poisson estimation using Goodman etc    ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

xSVD = [axhat(:); bxhat(:); kthat(:)];

NbParam=length(xSVD);

fprintf('Compute Poisson LL using Goodman. M: log(q) = a_x+b_x k_t\n')
x = xSVD;
dist=10^10;iterCtr=0;
while dist>10^(-10) && iterCtr<10000
    iterCtr=iterCtr+1;
    [x1,axhat1,bxhat1,kthat1] = LCGoodmanPoiss(x,DxtTrue,ExtTrue,nbAges,nbYrEst);
    % as described in Brouhns, Denuit, Vermunt
    dist = sum(abs(x1-x));
    x=x1;
end


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' get ready for the Renshaw-Haberman estimation    ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% number of cohorts:
nbCohorts = nbAges-Ns1+nbYrEst-Ns2-1;

% define a structure with general parameters
EstS.model = 'log(qxt)=a_x + b_x k_t + b2_x \gam_z, Estimate with Poisson ';
EstS.alpha = alpha;
EstS.omega = omega;
EstS.yearStart = yearStart;
EstS.yearEnd = yearEnd;
EstS.Ns1 = Ns1;
EstS.Ns2 = Ns2;
EstS.nbCohorts = nbCohorts;
EstS.nbYrEst = nbYrEst;
EstS.nbAges  = nbAges;


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' implement now the complete Poisson estimation with bilinear cohort effect   ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

gzhat=0.1*ones(nbCohorts,1); %some starting value if fmincon does not exist


bxhat2 = 0.01*ones(nbAges,1);
x0 = [axhat1(:);bxhat1(:);kthat1(:);bxhat2(:);gzhat(:)];

NbParam=length(x0);

tic
dist=10^10;iterCtr=0;
while dist>0.00001 && iterCtr<10000
    
    iterCtr=iterCtr+1;
    
    x1 = RH_GoodmanPoiss(x0,DxtTrue,ExtTrue,nbAges,nbCohorts,Ns1,Ns2,nbYrEst);
    
    %writedist(x0,x1,nbAges,nbcohorts,Ns1,Ns2,nbYrEst);
    
    % as described in Brouhns, Denuit, Vermunt
    dist = sum(abs(x1-x0));
    x0=x1;
end
EstS.tEst=toc;
fprintf('RH took %8.4f seconds\n',EstS.tEst)
fprintf('RH took %8.0f iterations\n',iterCtr)
EstS.dist=dist;
EstS.iterCtr=iterCtr;
xRH=x1;

plotRHestimates(xRH,nbAges,nbCohorts,Ns1,Ns2,nbYrEst);


%figure()
%plotResiduals(resid,yearv,age1,Ns1,Ns2);

%figure()
[qxMhat,Dxthat] = getRHqmatrix(xRH,ExtTrue,nbAges,nbCohorts,Ns1,Ns2,nbYrEst);

%plotqandD(qxM,qxMhat,agev);

resids=qxM-qxMhat;
EstS.resids=resids;

N  = nbAges; % to keep presentation nice below
T  = nbYrEst;
Nc = nbCohorts;

p1=1;    p2=N;     axhat1  = xRH(p1:p2);
p1=p2+1; p2=2*N;   bxhat1  = xRH(p1:p2);
p1=p2+1; p2=p2+T;  kthat1  = xRH(p1:p2);
p1=p2+1; p2=p2+N;  bxhat2  = xRH(p1:p2);
p1=p2+1; p2=p2+Nc; gzhat   = xRH(p1:p2);

EstS.axhat1 = axhat1;
EstS.bxhat1 = bxhat1;
EstS.kthat1 = kthat1;
EstS.bxhat2 = bxhat2;
EstS.gzhat  = gzhat;

qxt = getRHqx(EstS);
disp('qxt(1:5,1:5)=')
niceprint(qxt(1:5,1:5))
