% mainCreateStruct.m
% reads the Human Mortality Database for one country and saves the data in a
% Structure which gets stored with Matlab format
clc, close all, clear all

[year,age,mxt,qxt,lxt,dxt,ext] = loadData();


LifeTableS.year = year;
LifeTableS.age  = age;
LifeTableS.mxt  = mxt;
LifeTableS.qxt  = qxt;
LifeTableS.lxt  = lxt;
LifeTableS.dxt  = dxt;
LifeTableS.ext  = ext;

LifeOut = 'LifeTableS.mat' ;

save(LifeOut,'LifeTableS')