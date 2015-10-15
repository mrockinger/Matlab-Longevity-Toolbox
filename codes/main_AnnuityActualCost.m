% main_AnnuityActualCost
% scope: perform a comparison between life table approach to compute annuity
% cost and actual cost.
%
% computes values of all annuities for people retiring at age 65.
% assumes last individual of a cohort dies at age 95.
% does computations for as many years as possible.
% This allows investigating the consequence of changes in mortality on
% annuities
close all, clear all, clc, echo off

interestrategrid=[0,1,2,3,4,5]'/100; % the various interst rates to use
nbinterestrate=length(interestrategrid);

retage = 65; % Retirement age
omega  = 95; % maximal age

Infilename='LifeTableS.mat';
load(Infilename);

year = LifeTableS.year;
age  = LifeTableS.age;
qxt  = LifeTableS.qxt;

yearu   = unique(year); % sequence of years in database.
nbyears = length(yearu);

ageIs  = age(1:111); % the various ages 0 to 110+

qx     = reshape(qxt,size(qxt,1)/nbyears,nbyears);

nbcohorts=nbyears-(omega-retage);

% Annuity is computed with Periodic table approach, i.e. static mortality rates
retidx = find(ageIs==retage); % index from when on retirement starts
oldidx = find(ageIs==omega);  % index when last annuity gets paid

AM=zeros(nbcohorts,nbinterestrate);
for ictr=1:nbinterestrate
    for yrctr=1:nbcohorts
        
        interestrate=interestrategrid(ictr);
        
        % select mortality rate
        qxsel  = qx(retidx:oldidx,yrctr); %all mortalities from retirement on
        
        %%%%%%%%%%%%%%%%%%
        AM(yrctr,ictr) = annuity(qxsel,interestrate);
        %%%%%%%%%%%%%%%%%%
        
    end
end

% value of annuity for various years for individual age x.
% Annuity is computed with Periodic table approach, i.e. static mortality rates
AM2=zeros(nbcohorts,nbinterestrate);
nbelts=omega-retage+1; % number of elements in a given cohort
for ictr=1:nbinterestrate
    for yrctr=1:nbcohorts
        
        interestrate=interestrategrid(ictr);
        
        % select mortality rate
        qxsel=zeros(nbelts,1);
        % get mortality rates along the diagonal
        for k=1:nbelts
            qxsel(k)  = qx(retidx+k-1,yrctr+k-1); %all mortalities from retirement on
        end
        %%%%%%%%%%%%%%%%%%
        AM2(yrctr,ictr) = annuity(qxsel,interestrate);
        %%%%%%%%%%%%%%%%%%
        
    end
end

Symbol={'-k';'--k';':k';'-.k';'-k';'--k'};
Lwidth=[1;1;1;1;2;2];

%
h=figure(); iv=cell(nbinterestrate,1);
for ictr=1:nbinterestrate
    plot(yearu(1:nbcohorts),AM(:,ictr),Symbol{ictr},'Linewidth',Lwidth(ictr)); hold on
    iv{ictr}=['{\it i=}',num2str(interestrategrid(ictr))];
end


Symbol={'-r';'--r';':r';'-.r';'-r';'--r'};
Lwidth=[1;1;1;1;2;2];
for ictr=1:nbinterestrate
    plot(yearu(1:nbcohorts),AM2(:,ictr),Symbol{ictr},'Linewidth',Lwidth(ictr)); hold on
    iv{ictr}=['{\it i=}',num2str(interestrategrid(ictr))];
end
hold off
legend(iv,'Location','NorthWest')
titlestr=['Life Annuity Values comparison of methods, United States (bk=period,red=cohort)'];
title(titlestr)




% now generate a table in LaTex format

fprintf('&\\multicolumn{6}{c}{ %s }\\\\\n','US')
fprintf(' ')
for ictr=1:nbinterestrate
    fprintf('&%8.2f ',interestrategrid(ictr))
end
fprintf('\\\\\n ')

yearsel=yearu(1:nbcohorts);

for yrctr = 1:nbcohorts
    fprintf('%8.2f ',yearsel(yrctr));
    for ictr=1:nbinterestrate
        fprintf('&%8.2f ',AM(yrctr,ictr));
    end
    fprintf('\\\\\n ')
end

% now generate a table in LaTex format
fprintf('&\\multicolumn{6}{c}{ %s }\\\\\n','US')
fprintf(' ')
for ictr=1:nbinterestrate
    fprintf('&%8.2f ',interestrategrid(ictr))
end
fprintf('\\\\\n ')

yearsel=yearu(1:nbcohorts);

for yrctr = 1:nbcohorts
    fprintf('%8.2f ',yearsel(yrctr));
    for ictr=1:nbinterestrate
        fprintf('&%8.2f ',AM2(yrctr,ictr));
    end
    fprintf('\\\\\n ')
end
