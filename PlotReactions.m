function FigureObject=PlotReactions(CtrlVar,MUA,Reactions,FigureObject)


persistent FigReactions

if nargin==4
    FigReactions=FigureObject;
    FigReactions=figure(FigReactions) ; 
elseif isempty(FigReactions)
    FigReactions=figure; hold off
   
else
    try
        FigReactions=figure(FigReactions) ; hold off
    catch
        FigReactions=figure; hold off
    end
end


x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);


hold on


Two= numel(Reactions.ubvb)>0 && numel(Reactions.ubvb)>0;

if numel(Reactions.ubvb)>0
    
    if Two
        subplot(1,2,1) 
    end
    PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar) ; 
    QuiverColorGHG(x,y,Reactions.ubvb(1:MUA.Nnodes),Reactions.ubvb(MUA.Nnodes+1:end),CtrlVar);
    title('uv Reactions') ; 
end

if numel(Reactions.h)>0
    
    if Two
        subplot(1,2,2)
        PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar) ; 
        axis tight
    end
    
    [AreaTri]=TriAreaFE(MUA.coordinates,MUA.connectivity);
    Reactions.h= 3*Reactions.h*min(AreaTri)/max(abs(Reactions.h))/CtrlVar.PlotXYscale^2;
    
    I=Reactions.h>0; scatter(x(I)/CtrlVar.PlotXYscale,y(I)/CtrlVar.PlotXYscale,Reactions.h(I),'b','d','filled')
    I=Reactions.h<0; scatter(x(I)/CtrlVar.PlotXYscale,y(I)/CtrlVar.PlotXYscale,-Reactions.h(I),'r','c','filled')
    title('h Reactions') ; 
end


xlabel(CtrlVar.PlotsXaxisLabel) ; ylabel(CtrlVar.PlotsYaxisLabel); 

FigureObject=FigReactions;

end