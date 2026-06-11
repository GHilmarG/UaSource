
%%
[Lon,Lat,BasID] = ReadDrainageBasinsZwally();  
[xB,yB]=ll2xy(Lat,Lon);


hold on
for I=1:27
   Ind=I==BasID;
   plot(xB(Ind)/1000,yB(Ind)/1000,'k')
   text(mean(xB(Ind))/1000,mean(yB(Ind))/1000,num2str(I),'color','k') 
end

xlabel('xps (km)'); ylabel('yps (km)')
axis equal
hold off
axis off



