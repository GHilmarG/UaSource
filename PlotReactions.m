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
PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar) ; 

hold on

if numel(Reactions.ubvb)>0
    QuiverColorGHG(x,y,Reactions.ubvb(1:MUA.Nnodes),Reactions.ubvb(MUA.Nnodes+1:end),CtrlVar);
end

if numel(Reactions.h)>0
    
    [AreaTri]=TriAreaFE(MUA.coordinates,MUA.connectivity);
    Reactions.h= 3*Reactions.h*min(AreaTri)/max(abs(Reactions.h))/CtrlVar.PlotXYscale^2;
    
    I=Reactions.h>0; scatter(x(I)/CtrlVar.PlotXYscale,y(I)/CtrlVar.PlotXYscale,Reactions.h(I),'b','d','filled')
    I=Reactions.h<0; scatter(x(I)/CtrlVar.PlotXYscale,y(I)/CtrlVar.PlotXYscale,-Reactions.h(I),'r','c','filled')

end

title('Reactions') ; xlabel(CtrlVar.PlotsXaxisLabel) ; ylabel(CtrlVar.PlotsYaxisLabel); 

FigureObject=FigReactions;

end