function [year,age,mxt,qxt,lxt,dxt,ext] = loadData()
% loadData loads the Human Mortality Dataset (here both genders USA) and returns
% the various components thereoff. If you wish to use the full set of data,
% please register with the Human Mortality Project 

%% Processing LifeTables (female, male and total)


% Counting the number of lines
LifeTablefilename='bltper_1x1.txt';
fileID = fopen(LifeTablefilename,'r');
line1 = fgetl(fileID);
line2 = fgetl(fileID);
line3 = fgetl(fileID);

countLines=0;
while ~feof(fileID)
    fgetl(fileID);
    countLines=countLines+1;
end
fclose(fileID);

% Storing the data
year =zeros(countLines,1);
age = zeros(countLines,1);
mxt = zeros(countLines,1);
qxt = zeros(countLines,1);
axt = zeros(countLines,1);
lxt = zeros(countLines,1);
dxt = zeros(countLines,1);
Lxt = zeros(countLines,1);
Txt = zeros(countLines,1);
ext  =zeros(countLines,1);


fileID = fopen(LifeTablefilename,'r');
line1 = fgetl(fileID);
line2 = fgetl(fileID);
line3 = fgetl(fileID);
for t=1:countLines
    
    aline=fgetl(fileID);
    
    year(t) = str2num(aline(3:6));
    
    ageraw=str2num(aline(15:18));
    age(t) = ageraw;
    
    mxraw = str2num(aline(20:30));
    if isempty(mxraw); mxraw=NaN; end
    mxt(t) = mxraw;
    
    qxraw = str2num(aline(31:39));
    if isempty(qxraw); qxraw=NaN; end
    qxt(t) = qxraw;
    
    axraw = str2num(aline(40:45));
    if isempty(axraw); axraw=NaN; end
    axt(t) = axraw;
    
    lxraw = str2num(aline(46:53));
    if isempty(lxraw); lxraw=NaN; end
    lxt(t) = lxraw;
    
    dxraw = str2num(aline(54:61));
    if isempty(dxraw); dxraw=NaN; end
    dxt(t) = dxraw;
    
    Lxraw = str2num(aline(62:69));
    if isempty(Lxraw); Lxraw=NaN; end
    Lxt(t) = Lxraw;
    
    Txraw = str2num(aline(70:78));
    if isempty(Txraw); Txraw=NaN; end
    Txt(t) = Txraw;
    
    exraw = str2num(aline(79:85));
    if isempty(exraw); exraw=NaN; end
    ext(t) = exraw;
    
end
fclose(fileID);

end

