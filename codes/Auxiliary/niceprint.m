function niceprint(x)
[r,c]=size(x);
fmt = repmat('%20.6f',1,c);
fmt = strcat(fmt,' \n');
fprintf(fmt,x');