function ModifyColormap


% sets colormap to gray around 0

[t1,t2]=caxis ;
range=(t2-t1)*linspace(0,1,size(colormap,1))+t1 ;


temp=range ; temp(temp>0)=NaN;
[~,I1]=max(temp);
temp=range ; temp(temp<0)=NaN;
[~,I2]=min(temp);
I=[I1;I2];
temp=colormap;
temp(I,:)=temp(I,:)*0+0.9;
colormap(gca,temp)


end




