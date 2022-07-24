function ModifyColormap(GrayLevel,Ncol,options)

%
% resets colormap to gray for values at GrayLevel, and uses different colorscales for values below and above GrayLevel. 
%
% Examples:
%
%
%% 
%
%   [X,Y,Z] = peaks(500); figure ; contourf(X,Y,Z,20) ; colorbar ; ModifyColormap  ; 
%
%
%   figure ; [X,Y,Z] = peaks(500); contourf(X,Y,Z,20) ; colorbar ; ModifyColormap(ShowGrayLevel=false,Ncol=10);
%
%
%   sp221=subplot(2,2,1);
%   ModifyColormap(ShowGrayLevel=false,Ncol=2028,handle=sp221)  ; 
%


arguments
    GrayLevel (1,1) double = nan
    Ncol (1,1) double = nan 
    options.GrayLevel (1,1) double=0
    options.Ncol (1,1) double = 1024
    options.ShowGrayLevel logical=true
    options.handle double=[]
end

if isnan(Ncol)
    Ncol=options.Ncol;
end

if isnan(GrayLevel)
    GrayLevel=options.GrayLevel;
end


if isempty(options.handle)
    options.handle=gca ;
end


temp=colormap(options.handle,othercolor("YlOrRd9",Ncol)) ;

[t1,t2]=caxis ;
range=(t2-t1)*linspace(0,1,size(temp,1))+t1 ;



[~,iloc]=min(abs(range-GrayLevel));

if options.ShowGrayLevel
    N=2;
    I=iloc-N:iloc+N;
    I(I<1)=[];
    I(I>options.Ncol)=[];
else
    N=0;
end


PosColorscale="YlOrRd9";
NegColorscale="YlGnBu8";

temp(iloc+N:Ncol,:)=othercolor(PosColorscale,Ncol-(iloc+N-1)) ;
temp(1:iloc-N,:)=flipud(othercolor(NegColorscale,iloc-N)) ;

if options.ShowGrayLevel
    temp(I,:)=temp(I,:)*0+0.95;
end

colormap(options.handle,temp)


end




