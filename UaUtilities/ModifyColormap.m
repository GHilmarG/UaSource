function ModifyColormap(GrayLevel,Ncol)


% sets colormap to gray around 0
% and, optionally, a differnet colorscales above and below the zero level
%

if nargin==0
    GrayLevel=0;
end

if nargin<2
    Ncol=size(colormap,1); 
else
    colormap(othercolor("YlOrRd9",Ncol)) ; 
end


[t1,t2]=caxis ;

range=(t2-t1)*linspace(0,1,size(colormap,1))+t1 ;



[~,iloc]=min(abs(range-GrayLevel));

N=2;
I=[iloc-N:iloc+N];
I(I<1)=[];
temp=colormap;


PosColorscale="YlOrRd9";
NegColorscale="YlGnBu8";
temp(iloc+N:Ncol,:)=othercolor(PosColorscale,Ncol-(iloc+N-1)) ; 
temp(1:iloc-N,:)=flipud(othercolor(NegColorscale,iloc-N)) ; 

temp(I,:)=temp(I,:)*0+0.95;

colormap(gca,temp)


end




