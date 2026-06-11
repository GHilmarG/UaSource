function PlotReactions(CtrlVar,MUA,F,Reactions)

narginchk(4,4)

%
% To calculate Reactions, use CalculateReactions.m
%

x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);

Two= numel(Reactions.ubvb)>0 && numel(Reactions.h)>0;

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
    end
    PlotFEmesh(MUA.coordinates,MUA.connectivity,CtrlVar) ;
    axis tight


    [AreaTri]=TriAreaFE(MUA.coordinates,MUA.connectivity);
    Rh= 3*Reactions.h*min(AreaTri)/max(abs(Reactions.h))/CtrlVar.PlotXYscale^2;

    I=Rh>0; scatter(x(I)/CtrlVar.PlotXYscale,y(I)/CtrlVar.PlotXYscale,Rh(I),'b','d','filled')
    I=Rh<0; scatter(x(I)/CtrlVar.PlotXYscale,y(I)/CtrlVar.PlotXYscale,-Rh(I),'r','c','filled')
    title("$h$-reactions",Interpreter="latex") ;
    subtitle("negative reactions in red, positive in blue",interpreter="latex")


    if nargin> 3
        ah=Reactions.h./F.rho/F.dt;
        UaPlots(CtrlVar,MUA,F,ah,FigureTitle=" Reactions.h/(F.rho F.dt) Reactions ")
        title("$h$-reactions/($\rho \, \Delta t)$",Interpreter="latex") ;
        subtitle(sprintf("t=%g",F.time),Interpreter="latex") ;

    end

end


xlabel(CtrlVar.PlotsXaxisLabel) ; ylabel(CtrlVar.PlotsYaxisLabel);


end