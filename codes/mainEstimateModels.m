% main_EstimateModels.m
% runs all the estimations for the various models. These are LC1, LC2, RH, CBD1, CBD2 and CBD3.
close all, clear all, clc


retage=65; % Retirement age

Infilename='LifeTableS.mat';
load(Infilename);

year = LifeTableS.year;
age  = LifeTableS.age;
qxt  = LifeTableS.qxt;

yearu   = unique(year); % sequence of years in database.
nbyears = length(yearu);

ageIs  = age(1:111); % the various ages 0 to 110+ 

qx     = reshape(qxt,size(qxt,1)/nbyears,nbyears);

nbyear = length(yearu);

alpha=0;
omega=110;
year=1961; 

% % models for the full age range. Estimates parameters for just one year  
disp('Heligman-Pollard Parametric model 1')
EstS =  HPmodel1(year,yearu,alpha,omega,qx); 
fprintf('%s\n',EstS.model);
fprintf('A=%10.4f\n',EstS.paramv(1));
fprintf('B=%10.4f\n',EstS.paramv(2));
fprintf('C=%10.4f\n',EstS.paramv(3));
fprintf('D=%10.4f\n',EstS.paramv(4));
fprintf('E=%10.4f\n',EstS.paramv(5));
fprintf('F=%10.4f\n',EstS.paramv(6));
fprintf('G=%10.4f\n',EstS.paramv(7));
fprintf('H=%10.4f\n',EstS.paramv(8));
fprintf('Exitflag%3.0f\n',EstS.exitflag);
fprintf('Root Mean Squared Errors %10.5f\n', RMSE(EstS.resids))
fprintf('Rsquare...................... %10.5f\n', EstS.Rsquare)


disp('Heligman-Pollard Parametric model 1.a')
EstS =  HPmodel1a(year,yearu,alpha,omega,qx); 
fprintf('%s\n',EstS.model);
fprintf('A=%10.4f\n',EstS.paramv(1));
fprintf('B=%10.4f\n',EstS.paramv(2));
fprintf('C=%10.4f\n',EstS.paramv(3));
fprintf('D=%10.4f\n',EstS.paramv(4));
fprintf('E=%10.4f\n',EstS.paramv(5));
fprintf('F=%10.4f\n',EstS.paramv(6));
fprintf('G=%10.4f\n',EstS.paramv(7));
fprintf('H=%10.4f\n',EstS.paramv(8));
fprintf('Exitflag%3.0f\n',EstS.exitflag);
fprintf('Root Mean Squared Errors %10.5f\n', RMSE(EstS.resids))
fprintf('Rsquare...................... %10.5f\n', EstS.Rsquare)


disp('Heligman-Pollard Parametric model 2')
EstS =  HPmodel2(year,yearu,alpha,omega,qx); 
fprintf('%s\n',EstS.model);
fprintf('A=%10.4f\n',EstS.paramv(1));
fprintf('B=%10.4f\n',EstS.paramv(2));
fprintf('C=%10.4f\n',EstS.paramv(3));
fprintf('D=%10.4f\n',EstS.paramv(4));
fprintf('E=%10.4f\n',EstS.paramv(5));
fprintf('F=%10.4f\n',EstS.paramv(6));
fprintf('G=%10.4f\n',EstS.paramv(7));
fprintf('H=%10.4f\n',EstS.paramv(8));
fprintf('K=%10.4f\n',EstS.paramv(9));
fprintf('Exitflag%3.0f\n',EstS.exitflag);
fprintf('Root Mean Squared Errors %10.5f\n', RMSE(EstS.resids))
fprintf('Rsquare...................... %10.5f\n', EstS.Rsquare)

disp('Heligman-Pollard Parametric model 3')
EstS =  HPmodel3(year,yearu,alpha,omega,qx); 
fprintf('%s\n',EstS.model);
fprintf('A=%10.4f\n',EstS.paramv(1));
fprintf('B=%10.4f\n',EstS.paramv(2));
fprintf('C=%10.4f\n',EstS.paramv(3));
fprintf('D=%10.4f\n',EstS.paramv(4));
fprintf('E=%10.4f\n',EstS.paramv(5));
fprintf('F=%10.4f\n',EstS.paramv(6));
fprintf('G=%10.4f\n',EstS.paramv(7));
fprintf('H=%10.4f\n',EstS.paramv(8));
fprintf('K=%10.4f\n',EstS.paramv(9));
fprintf('Exitflag%3.0f\n',EstS.exitflag);
fprintf('Root Mean Squared Errors %10.5f\n', RMSE(EstS.resids))
fprintf('Rsquare...................... %10.5f\n', EstS.Rsquare)

% % models that take the evolution over time somehow into account
yearStart = 1947;
yearEnd   = 2010;
alpha     = 55;
omega     = 85;

disp('Lee and Carter under Gaussian innovations')
EstS =  LCmodel1(yearStart,yearEnd,yearu,alpha,omega,ageIs,qx);
fprintf('Root Mean Squared Errors %10.5f\n', RMSE(EstS.resids))


disp('Lee and Carter under Poisson innovations')
EstS =  LCmodel2(yearStart,yearEnd,yearu,alpha,omega,ageIs,qx);
fprintf('%s\n',EstS.model);
fprintf('Required time for ML optimization %10.5f\n', EstS.tEst)
fprintf('Final distance between parameter vectors %12.11f\n', EstS.dist)
fprintf('Required number of iterations %10.5f\n', EstS.iterCtr)
fprintf('Root Mean Squared Errors %10.5f\n', RMSE(EstS.resids))

disp('Renshaw and Haberman under Poisson innovations')
EstS = RHmodel(yearStart,yearEnd,yearu,alpha,omega,ageIs,qx);
fprintf('%s\n',EstS.model);
fprintf('Required time for ML optimization %10.5f\n', EstS.tEst)
fprintf('Final distance between parameter vectors %12.11f\n', EstS.dist)
fprintf('Required number of iterations %10.5f\n', EstS.iterCtr)
fprintf('Root Mean Squared Errors %10.5f\n', RMSE(EstS.resids))


disp('Cairns, Blake and Dowd model 1 under Poisson innovations')
EstS =  CBDmodel1(yearStart,yearEnd,yearu,alpha,omega,ageIs,qx);
fprintf('%s\n',EstS.model);
fprintf('Required time for ML optimization %10.5f\n', EstS.tEst)
fprintf('Final distance between parameter vectors %12.11f\n', EstS.dist)
fprintf('Required number of iterations %10.5f\n', EstS.iterCtr)
fprintf('Root Mean Squared Errors %10.5f\n', RMSE(EstS.resids))


disp('Cairns, Blake and Dowd model 2 under Poisson innovations')
EstS =  CBDmodel2(yearStart,yearEnd,yearu,alpha,omega,ageIs,qx);
fprintf('%s\n',EstS.model);
fprintf('Required time for ML optimization %10.5f\n', EstS.tEst)
fprintf('Final distance between parameter vectors %12.11f\n', EstS.dist)
fprintf('Required number of iterations %10.5f\n', EstS.iterCtr)
fprintf('Root Mean Squared Errors %10.5f\n', RMSE(EstS.resids))


disp('Cairns, Blake and Dowd model 3 under Poisson innovations')
EstS =  CBDmodel3(yearStart,yearEnd,yearu,alpha,omega,ageIs,qx);
fprintf('%s\n',EstS.model);
fprintf('Required time for ML optimization %10.5f\n', EstS.tEst)
fprintf('Final distance between parameter vectors %12.11f\n', EstS.dist)
fprintf('Required number of iterations %10.5f\n', EstS.iterCtr)
fprintf('Root Mean Squared Errors %10.5f\n', RMSE(EstS.resids))