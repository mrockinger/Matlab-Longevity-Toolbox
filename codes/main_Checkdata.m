% main_Checkdata
% loads the data and prints the data to a checkfile
% Name of checkfile is "CheckData.txt"

clc, close all, clear all

[year,age,mxt,qxt,lxt,dxt,ext] = loadData();


NbElts=size(year,1);

fOut=fopen('CheckData.txt','w');

for t=1:NbElts
   fprintf(fOut,'%5.0f %4.0f %8.5f %8.5f %7.0f %5.0f %6.2f \n',year(t),age(t),mxt(t),qxt(t),lxt(t),dxt(t),ext(t)); 
end

fclose(fOut);