function EstS =  HPmodel2(year,yearu,alpha,omega,qx)

% USAGE :
% Estimates the parameters of the 9 parameter Heligman Pollard model

% INPUTS :
% year the year for which want to obtain the estimates
% yearu the year corresponsing to the mortality matrix
% alpha, omega are 0 and 110 for this model since cover entire age range
% qx the mortality rates
%
%OUTPUTS :
% EstS.model - name of model estimated
% EstS.resids - difference between mortality and estimates
% EstS.paramv - vector of parameters A,B,... as in note
% EstS.exitflag - a flag to indicate the convergence criterion of optimization.
% EstS.resids - residuals defined as the mortality minus estimated mortality
% EstS.Rsquare - a measure of goodness of fit

EstS.model='HPmodel2';

NbParam=9; % this is the 9 parameter model

yidx = find(year==yearu); % year for which one wishes to perform the parameter estimation
ages = (alpha:omega)'; % full age range

yrs  = yearu(yidx); % verification
qxts = qx(:,yidx);

ages = ages(1:end-1); % remove 110 year olds
qxts = qxts(1:end-1);

figure(1)
plot(ages,log(qxts),'ok'); hold on
title('mortality rates')

%estimate Heligman and Pollard. Model 1a. Parameter Nbr=8.
scal = ones(NbParam,1);
scal(1) = 1/1000;
scal(2) = 1/100;
scal(3) = 1/10;
scal(4) = 1/10000;
scal(5) = 10;
scal(6) = 10;
scal(7) = 1/100000;
scal(8) = 1;
scal(9) = 1;

%estimate Heligman and Pollard. Model 2. Parameter Nbr=9. senescense
%correction in denuminator


A = [];
b = [];
Aeq = [];
beq = [];
lb = 0.000001*ones(NbParam,1);
ub = [500;1000;10;150;2.5;10;80; 1.3;3];

nonlcon = [];

%Level of display. 
%'off' displays no output; 
%'iter' displays output at each iteration; 
%'final' displays just the final output; 
%'notify' (default) displays output only if the function does not converge. 
options = optimset('display','final','MaxFunEvals',10000,...
                   'TolCon',0.0000000001,...
                   'TolFun',0.0000000001,...
                   'TolX',0.000000001,...
                   'MaxIter',10000);
x01=[ 7.7983; 2.0776; 1.7867; 16.7963; 0.8303; 2.2038; 19.1715; 1.0804;1];



%  [x,fval,exitflag,output] = fminsearchbnd(fun,x0,LB,UB,options,varargin)
[paramv,fval,exitflag,output] = fminsearchbnd(@HP2,x01,lb,ub,options,ages,qxts,scal);
EstS.paramv   = scal.*paramv;
EstS.exitflag = exitflag;

y=HP2_Shape(paramv,ages,scal);

% construct an Rsquare
yi=log(qxts);
fi=log(y);
resids = yi-fi;
EstS.resids=resids;


SStot = sum( (yi-mean(yi)).^2 );
SSres = sum( (yi-fi).^2 );
EstS.Rsquare = 1-SSres/SStot;
fprintf('Rsquare is %8.4f\n', EstS.Rsquare)

plot(ages,log(y),'-k'); hold off

disp('True parameters')
alphabet=['A','B','C','D','E','F','G','H','K'];
for paridx=1:NbParam
  str=[alphabet(paridx),'= %8.4f\n'];
  fprintf(str,scal(paridx)*paramv(paridx))
end
