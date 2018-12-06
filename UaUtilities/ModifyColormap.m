function ModifyColormap


% sets colormap to gray around 0

[t1,t2]=caxis ;
range=(t2-t1)*linspace(0,1,size(colormap,1))+t1 ;
%I=abs(range)<0.3;
[~,I]=min(abs(range));
temp=colormap;
temp(I,:)=temp(I,:)*0+0.9;
colormap(gca,temp)


end




