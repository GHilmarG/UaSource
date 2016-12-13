function AdjointResultsPlots(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo)

save TestSaveAdjointResults
%%
CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
%load TestSave
us=F.ub+F.ud; vs=F.vb+F.vd;

x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);

GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
xGL=[] ; yGL=[] ;

%%
figure

subplot(2,2,1)

PlotMeshScalarVariable(CtrlVar,MUA,Meas.us) ; hold on ;
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
title('us Meas on numerical grid') ;

subplot(2,2,2)

PlotMeshScalarVariable(CtrlVar,MUA,Meas.vs) ; hold on ;
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
title('vs Meas on numerical grid') ;

usError=sqrt(spdiags(Meas.usCov));
vsError=sqrt(spdiags(Meas.vsCov));
wsError=sqrt(spdiags(Meas.wsCov));

subplot(2,2,3)


PlotMeshScalarVariable(CtrlVar,MUA,usError) ; hold on ;
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
title('us error on numerical grid') ;

subplot(2,2,4)
PlotMeshScalarVariable(CtrlVar,MUA,vsError) ; hold on ;
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
title('vs error on numerical grid') ;



%%
if strcmp(CtrlVar.AdjointGrad,'A')
    
    figure ; PlotMeshScalarVariable(CtrlVar,MUA,log10(InvFinalValues.AGlen));
    title('log10(InvFinalValues.AGlen)') ; cbar=colorbar; title(cbar, '(a^{-1} kPa^{-3})');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    
    figure ; PlotMeshScalarVariable(CtrlVar,MUA,log10(InvStartValues.AGlen));
    title('log10(Astart)') ; cbar=colorbar; title(cbar, '(a^{-1} kPa^{-3})');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    
    figure ; PlotMeshScalarVariable(CtrlVar,MUA,log10(InvFinalValues.AGlen)-log10(InvStartValues.AGlen));
    title('log10(InvFinalValues.AGlen)-log10(InvStartValues.AGlen)') ; cbar=colorbar; title(cbar, '(a^{-1} kPa^{-3})');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
end

%%
if strcmp(CtrlVar.AdjointGrad,'C')
    
    figure ; PlotMeshScalarVariable(CtrlVar,MUA,log10(InvFinalValues.C));
    title('log10(InvFinalValues.C)') ; cbar=colorbar; title(cbar, '(m/a/kPa^m)');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    
    figure ; PlotMeshScalarVariable(CtrlVar,MUA,log10(InvStartValues.C));
    title('log10(Cstart)') ; cbar=colorbar; title(cbar, '(m/a/kPa^m)');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    
    figure ; PlotMeshScalarVariable(CtrlVar,MUA,log10(InvFinalValues.C)-log10(InvStartValues.C));
    title('log10(InvFinalValues.C)-log10(Cstart)') ; cbar=colorbar; title(cbar, '(m/a/kPa^m)');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
end


%%
[tbx,tby,tb,beta2] = CalcBasalTraction(CtrlVar,MUA,F.ub,F.vb,InvFinalValues.C,InvFinalValues.m,GF) ;
figure
PlotMeshScalarVariable(CtrlVar,MUA,tb) ;
title(' tb ') ; cbar=colorbar; title(cbar, '(kPa)');

%%
CtrlVar.VelPlotIntervalSpacing='log10';
figure
subplot(2,2,1);
QuiverColorGHG(x,y,(us-Meas.us)./usError,(vs-Meas.vs)./vsError,CtrlVar);
title('((us-Meas.us)/usError,(vs-Meas.vs)/vsError)') ;
hold on
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)

subplot(2,2,2);
QuiverColorGHG(x,y,us-Meas.us,vs-Meas.vs,CtrlVar); axis equal ; title('(us-Meas.us,v-Meas.vs)') ;
hold on ;
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)

subplot(2,2,3);
QuiverColorGHG(x,y,Meas.us,Meas.vs,CtrlVar); axis equal ;
title('(Meas.us,Meas.vs)') ;
hold on ;
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)

subplot(2,2,4);
QuiverColorGHG(x,y,us,vs,CtrlVar); axis equal ; title('(us,vs)') ;
hold on ;
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)

%%

IFigGradients=figure('Name','Gradients','NumberTitle','off');
subplot(3,1,1) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdp) ; title('dJdp')
subplot(3,1,2) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dIdp) ; title('dIdp')
subplot(3,1,3) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dRdp) ; title('dRdp')
%IFigGradients.Position=[1098.7 638.71 1096 518.29];
%%
CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
if ~isempty(Priors.TrueC)
    tFig1=figure('Name','True and estimated','NumberTitle','off');
    subplot(1,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,Priors.TrueC) ; title('True C')
    hold on ; PlotMuaMesh(CtrlVar,MUA,[],'w')
    subplot(1,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.C) ; title('Retrieved C')
    hold on ; PlotMuaMesh(CtrlVar,MUA,[],'w')
    tFig1.Units='normalized';
    tFig1.Position=[0.5 0.51 0.5 0.4];
    
    tFig2=figure('Name','Difference between true and estimated','NumberTitle','off');
    %PlotMeshScalarVariable(CtrlVar,MUA,Priors.TrueC-InvFinalValues.C);
    
    subplot(1,3,1)
    tri=PlotNodalVariableAsTriSurface(CtrlVar,MUA,[],InvFinalValues.C);
    xlabel('x') ; ylabel('y') ; title('Inverted slipperiness')
    colorbar('south')
    
    subplot(1,3,2)
    PlotNodalVariableAsTriSurface(CtrlVar,MUA,tri,Priors.TrueC);
    xlabel('x') ; ylabel('y') ; title('True slipperiness')
    colorbar('south')
    
    subplot(1,3,3)
    PlotNodalVariableAsTriSurface(CtrlVar,MUA,tri,InvFinalValues.C-Priors.TrueC);
    xlabel('x') ; ylabel('y') ; title('True slipperiness')
    title('Slipperiness: True-Estimated')
    colorbar('south')

    tFig2.Units='normalized';
    tFig2.Position=[0.1 0.2 0.8 0.5];
end

%%
if ~isempty(RunInfo.Inverse.J)
    figure
    yyaxis left
    semilogy(RunInfo.Inverse.Iterations,RunInfo.Inverse.J,'-bo','LineWidth',2)
    ylabel('J')
    
    if ~isempty(RunInfo.Inverse.I) && ~isempty(RunInfo.Inverse.R)
        
        hold on
        semilogy(RunInfo.Inverse.Iterations,RunInfo.Inverse.I,'-gx')
        ylabel('J & I')
        hold off
        yyaxis right
        semilogy(RunInfo.Inverse.Iterations,RunInfo.Inverse.R,'-r+')
        ylabel('R')
        legend('Objective function','I','R','Location','southwest')
    else
        legend('Objective function')
    end
    
    yyaxis left
    xlabel('Inverse iteration') ;
    hold off
end
%%
end
