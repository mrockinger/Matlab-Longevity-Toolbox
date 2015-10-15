% main_AnnuityCost
% computes the cost of a life annuity over the various years and various
% interest rates. Plots those series over time and make a table with the costs
% in LaTex format. 

close all, clear all, clc

interestrategrid=[0,1,2,3,4,5]'/100; % the various interst rates to use
nbinterestrate=length(interestrategrid);

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

%value of annuity for various years for individual age x.
% Annuity is computed with Periodic table approach, i.e. static mortality rates
AM=zeros(nbyear,nbinterestrate);
for ictr=1:nbinterestrate
    for yrctr=1:nbyear
        
        interestrate=interestrategrid(ictr);
        
        % select mortality rate        
        retidx = find(ageIs==retage); % index from when on retirement starts
        qxsel  = qx(retidx:end,yrctr); %all mortalities from retirement on
        
        %%%%%%%%%%%%%%%%%%
        AM(yrctr,ictr) = annuity(qxsel,interestrate);
        %%%%%%%%%%%%%%%%%%
        
    end
end

% formats for the plot
Symbol={'-k';'--k';':k';'-.k';'-b';'--b';':b';'-.b';'-k';'--k';':k';'-.k';'-b';'--b';':b';'-.b'};
Lwidth=[1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2;3;3;3;3];


h=figure(); iv=cell(nbinterestrate,1);
for ictr=1:nbinterestrate
    plot(yearu,AM(:,ictr),Symbol{ictr},'Linewidth',Lwidth(ictr)); hold on
    iv{ictr}=['{\it i=}',num2str(interestrategrid(ictr))];
end
legend(iv,'Location','NorthWest')
titlestr='Life Annuity Values (period life table approach), United States ';
title(titlestr)



