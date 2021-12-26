function ModifyColormap


% sets colormap to gray around 0
% and, optionally, a differnet colorscales above and below the zero level
%

[t1,t2]=caxis ;
Ncol=size(colormap,1); 
range=(t2-t1)*linspace(0,1,size(colormap,1))+t1 ;


[~,iloc]=min(abs(range));

N=2;
I=[iloc-N:iloc+N];
temp=colormap;


PosColorscale="YlOrRd9";
NegColorscale="YlGnBu8";
temp(iloc+N:Ncol,:)=othercolor(PosColorscale,Ncol-(iloc+N-1)) ; 
temp(1:iloc-N,:)=flipud(othercolor(NegColorscale,iloc-N)) ; 

temp(I,:)=temp(I,:)*0+0.95;

colormap(gca,temp)


end




