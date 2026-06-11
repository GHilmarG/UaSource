%%

function [xps,yps,lat,log]=ReadDepoorterGroundingLine

M=dlmread('Antarctic_grounding_line.tab','\t',14,0);
lat=M(:,1);  lon=M(:,2);  ID=M(:,3) ;
 
[xps,yps]=ll2xy(lat,lon);


%figure ; plot(xps,yps,'-') ; axis equal


end