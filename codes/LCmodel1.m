function EstS =  LCmodel1(yearStart,yearEnd,yearv,alpha,omega,agev,qx)


% USAGE : 
% Estimates the parameters of the LC model using SVD with no Poisson
% adjustment 

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
[nbAges,nbYrEst] = size(qxM);

% If there is a zero in the matrix can lead sometimes to problems
e=qxM<0.0000001; qxM(e)=0.00001;

% define a structure with general parameters
EstS.model = 'log(qxt)=a_x + b_x k_t Estimate with SVD ';
EstS.alpha = alpha;
EstS.omega = omega;
EstS.yearStart = yearStart;
EstS.yearEnd = yearEnd;
EstS.nbYrEst = nbYrEst;
EstS.nbAges  = nbAges;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' get the Lee-Carter 1992 decomposition')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
tic
[axhat,bxhat,kthat,qxMhat] = LC1(qxM);

EstS.tEst=toc;

fprintf('SVD took %8.4f seconds\n',EstS.tEst)

EstS.axhat=axhat;
EstS.bxhat=bxhat;
EstS.kthat=kthat;

EstS.resids = qxM-qxMhat;

agec   = agev(retidx:oldidx); % age and year conform with dimension of q matrix
yearc  = yearv(ySidx:yEidx);  % 


%figure()
%plotResiduals(resids,yearc,agec,10,10);

printParameters(axhat,bxhat,kthat,agec,yearc,qxM);
    
figure()
plotqandD(qxM,qxMhat,agec);




