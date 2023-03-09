function cmap=ModifyColormap(GrayLevel,Ncol,options)

%
% resets colormap to gray for values at GrayLevel, and uses different colorscales for values below and above GrayLevel. 
%
%
%%

%%
% Examples:
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
% Modify existing colormap so that values over the range of +/-100 around 0 are set to gray
%
%   ModifyColormap(ChangeColormap=false,GrayLevel=0,GrayLevelRange=100);
%
% Use with cmocean, setting colormap to gray around -100 to 100
%
%   CM=cmocean('balanced',25,'pivot',0) ; colormap(CM); ModifyColormap(100,nan,ChangeColormap=false) ;
%
%%

arguments
    GrayLevel (1,1) double = nan
    Ncol (1,1) double = nan 
    options.GrayLevel (1,1) double=0
    options.GrayLevelRange (1,1) double=nan
    options.Ncol (1,1) double = 1024
    options.ShowGrayLevel logical=true
    options.handle double=[]
    options.ChangeColormap=true; 
end

if isnan(Ncol)
    Ncol=options.Ncol;
end

if isnan(GrayLevel)
    GrayLevel=options.GrayLevel;
else
    options.GrayLevelRange=GrayLevel;
end


if isempty(options.handle)
    options.handle=gca ;
end

if options.ChangeColormap
    cmap=colormap(options.handle,othercolor("YlOrRd9",Ncol)) ;
else
    cmap=colormap;
end

Ncol=size(cmap,1); 

[t1,t2]=clim ;
range=(t2-t1)*linspace(0,1,size(cmap,1))+t1 ;



[~,iloc]=min(abs(range-GrayLevel));

if options.ShowGrayLevel  % set colors to gray around GrayLevel value over the range of values 
                          % as specified by GrayLevelRangeq
    if isnan(options.GrayLevelRange)
        
        N=round(size(cmap,1)*0.05) ;  % if range not specified, to 5% of total range 

    else
        N=round(options.GrayLevelRange/(range(2)-range(1)));
    end
    I=iloc-N:iloc+N;
    I(I<1)=[];
    I(I>Ncol)=[];
else
    N=0;
end

if options.ChangeColormap  % set different colormaps for pos and neg values
    PosColorscale="YlGnBu8";
    NegColorscale="YlOrRd9";
    

    cmap(iloc+N:Ncol,:)=othercolor(PosColorscale,Ncol-(iloc+N-1)) ;
    cmap(1:iloc-N,:)=flipud(othercolor(NegColorscale,iloc-N)) ;
end


if options.ShowGrayLevel  % set colors around zero values to gray
    cmap(I,:)=cmap(I,:)*0+0.95;
end

colormap(options.handle,cmap)


end




