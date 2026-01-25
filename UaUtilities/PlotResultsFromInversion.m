




function PlotResultsFromInversion(UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo)



%%
% function PlotResultsFromInversion(UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo)
%
%
% Does what it says on the tin.
%
%  Example:
%
% load InversionRestartFile
% PlotResultsFromInversion(UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo)
%
%
% It is also possible to enter the name of the restart file as the first, and only, argument. Then the restart file will be
% first loaded, and then plotted.
%
%
% Note: This function is used by Ua for plotting results from an inversion.
%
%%

if isstring(UserVar) && (isfile(UserVar)  || isfile(UserVar+".mat"))

    fprintf("loading and plotting results from %s \n",UserVar)

    load(UserVar,"UserVarInRestartFile","CtrlVarInRestartFile","MUA","BCs","F","BCsAdjoint","InvStartValues","InvFinalValues","Priors","Meas","RunInfo") ;

    CtrlVar=CtrlVarInRestartFile;
    UserVar=UserVarInRestartFile;

end


%%


CtrlVar.MUA.MassMatrix=true;
MUA=UpdateMUA(CtrlVar,MUA);


%%

fprintf(' Plotting results from inversion...')

CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
us=F.ub+F.ud; vs=F.vb+F.vd;

x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);

GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,F.GF,CtrlVar); xGL=[] ; yGL=[] ;

%%




figMeas=FindOrCreateFigure('Measurements') ; clf(figMeas)

T=tiledlayout("flow");

nexttile

cbar=UaPlots(CtrlVar,MUA,F,Meas.us,CreateNewFigure=false);

xlabel(CtrlVar.PlotsXaxisLabel,'interpreter','latex');
ylabel(CtrlVar.PlotsYaxisLabel,'interpreter','latex');
title('us Meas on numerical grid') ;
subtitle("")

nexttile

cbar=UaPlots(CtrlVar,MUA,F,Meas.vs,CreateNewFigure=false);
xlabel(CtrlVar.PlotsXaxisLabel,'interpreter','latex');
ylabel(CtrlVar.PlotsYaxisLabel,'interpreter','latex');
title('vs Meas on numerical grid') ;
subtitle("")

if ~isempty(Meas.dhdt)  && contains(CtrlVar.Inverse.Measurements,"-dhdt")

    nexttile
    cbar=UaPlots(CtrlVar,MUA,F,Meas.dhdt,CreateNewFigure=false);

    xlabel(CtrlVar.PlotsXaxisLabel,'interpreter','latex');
    ylabel(CtrlVar.PlotsYaxisLabel,'interpreter','latex');
    title('dh/dt Meas on numerical grid') ;
    subtitle("")
end

usError=sqrt(spdiags(Meas.usCov));
vsError=sqrt(spdiags(Meas.vsCov));
dhdtError=sqrt(spdiags(Meas.dhdtCov));



nexttile
cbar=UaPlots(CtrlVar,MUA,F,usError,CreateNewFigure=false);
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
title('us error on numerical grid') ;
subtitle("")

nexttile
cbar=UaPlots(CtrlVar,MUA,F,vsError,CreateNewFigure=false);
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
title('vs error on numerical grid') ;
subtitle("")

if ~isempty(Meas.dhdt)  && contains(CtrlVar.Inverse.Measurements,"-dhdt")

    nexttile
    cbar=UaPlots(CtrlVar,MUA,F,dhdtError,CreateNewFigure=false);
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    title('dh/dt error on numerical grid') ;
    subtitle("")

end
T.Padding="tight";   T.TileSpacing="tight";


%%
if contains(upper(CtrlVar.Inverse.InvertFor),'A')

    figAeI=FindOrCreateFigure('A at the end of inversion') ; clf(figAeI)
    PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.AGlen);
    set(gca,'ColorScale','log')
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    CtrlVar.PlotNodes=0 ; % PlotMuaMesh(CtrlVar,MUA,[],'k') ;
    title("$A$ at end of inversion",Interpreter="latex")
    cbar=colorbar; title(cbar, '($\mathrm{a}^{-1}$ $\mathrm{kPa}^{-3}$)',interpreter="latex");
    colormap(othercolor("Mtemperaturemap",1028))
    PlotMuaBoundary(CtrlVar,MUA,'k');
    ColorbarLimits=10.^[mean(log10(InvFinalValues.AGlen))-4*std(log10(InvFinalValues.AGlen))  mean(log10(InvFinalValues.AGlen))+4*std(log10(InvFinalValues.AGlen))];
    if ColorbarLimits(1)==ColorbarLimits(2)
        Eps=10*eps(ColorbarLimits(1));
        ColorbarLimits(1)=ColorbarLimits(1)-Eps;
        ColorbarLimits(2)=ColorbarLimits(2)+Eps;
    end
    clim(ColorbarLimits)

    figAsI=FindOrCreateFigure('A at the start of inversion') ; clf(figAsI)
    PlotMeshScalarVariable(CtrlVar,MUA,InvStartValues.AGlen);
    set(gca,'ColorScale','log')
    title("$A$ at start of inversion",Interpreter="latex")
    cbar=colorbar; title(cbar, '($\mathrm{a}^{-1}$ $\mathrm{kPa}^{-3}$)',interpreter="latex");
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    colormap(othercolor("Mtemperaturemap",1028))
    PlotMuaBoundary(CtrlVar,MUA,'k');
    clim(ColorbarLimits)

    figAcI=FindOrCreateFigure('Change in A during inversion run') ; clf(figAcI)
    PlotMeshScalarVariable(CtrlVar,MUA,log10(InvFinalValues.AGlen)-log10(InvStartValues.AGlen));
    title('log10(InvFinalValues.AGlen)-log10(InvStartValues.AGlen)') ;
    cbar=colorbar; title(cbar, '($\mathrm{a}^{-1}$ $\mathrm{kPa}^{-3}$)',interpreter="latex");
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel,Interpreter="latex")  ; ylabel(CtrlVar.PlotsYaxisLabel,Interpreter="latex")
    CM=cmocean('balanced',25,'pivot',0) ; colormap(CM);
    PlotMuaBoundary(CtrlVar,MUA,'k');


end

%%
if contains(upper(CtrlVar.Inverse.InvertFor),'C')

    ColorbarLimits=10.^[mean(log10(InvFinalValues.C))-4*std(log10(InvFinalValues.C))  mean(log10(InvFinalValues.C))+4*std(log10(InvFinalValues.C))];
    if ColorbarLimits(1)==ColorbarLimits(2)
        Eps=10*eps(ColorbarLimits(1));
        ColorbarLimits(1)=ColorbarLimits(1)-Eps;
        ColorbarLimits(2)=ColorbarLimits(2)+Eps;
    end
    figCeI=FindOrCreateFigure('C at the end of inversion') ; clf(figCeI)
    PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.C);
    set(gca,'ColorScale','log')
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    clim(ColorbarLimits)

    xlabel(CtrlVar.PlotsXaxisLabel,Interpreter="latex")  ; ylabel(CtrlVar.PlotsYaxisLabel,Interpreter="latex")
    CtrlVar.PlotNodes=0 ; % PlotMuaMesh(CtrlVar,MUA,[],'k') ;
    title("$C$ at end of inversion",Interpreter="latex")
    cbar=colorbar; title(cbar, '($\mathrm{m}\,\mathrm{yr}^{-1}\,\mathrm{kPa}^{-m}$)','interpreter','latex');
    colormap(othercolor("Mtemperaturemap",1028))
    PlotMuaBoundary(CtrlVar,MUA,'k');
    clim(ColorbarLimits)

    figCbI=FindOrCreateFigure('C at the beginning of inversion') ; clf(figCbI)
    PlotMeshScalarVariable(CtrlVar,MUA,InvStartValues.C);
    set(gca,'ColorScale','log')
    title("$C$ at start of inversion",Interpreter="latex")
    cbar=colorbar; title(cbar, '($\mathrm{m}\,\mathrm{yr}^{-1}\,\mathrm{kPa}^{-m}$)','interpreter','latex');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel,Interpreter="latex")  ; ylabel(CtrlVar.PlotsYaxisLabel,Interpreter="latex")
    colormap(othercolor("Mtemperaturemap",1028))
    PlotMuaBoundary(CtrlVar,MUA,'k');
    clim(ColorbarLimits)

    figCcI=FindOrCreateFigure('Change in C during inversion run') ; clf(figCcI)
    PlotMeshScalarVariable(CtrlVar,MUA,log10(InvFinalValues.C)-log10(InvStartValues.C));
    title('log10(InvFinalValues.C)-log10(Cstart)') ;
    cbar=colorbar; title(cbar, '($\mathrm{m}\,\mathrm{yr}^{-1}\,\mathrm{kPa}^{-m}$)','interpreter','latex');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel,Interpreter="latex")  ; ylabel(CtrlVar.PlotsYaxisLabel,Interpreter="latex")
    CM=cmocean('balanced',25,'pivot',0) ; colormap(CM);
    PlotMuaBoundary(CtrlVar,MUA,'k');
end

if contains(CtrlVar.Inverse.InvertFor,'b')

    figure ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.b);
    title('InvFinalValues.b') ; cbar=colorbar; title(cbar, '(m)');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel,Interpreter="latex")  ; ylabel(CtrlVar.PlotsYaxisLabel,Interpreter="latex")


    figure ; PlotMeshScalarVariable(CtrlVar,MUA,InvStartValues.b);
    title('bstart') ; cbar=colorbar; title(cbar, '(m)');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel,Interpreter="latex")  ; ylabel(CtrlVar.PlotsYaxisLabel,Interpreter="latex")


    figure ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.b-InvStartValues.b);
    title('InvFinalValues.b-bstart') ; cbar=colorbar; title(cbar, '(m)');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel,Interpreter="latex")  ; ylabel(CtrlVar.PlotsYaxisLabel,Interpreter="latex")

    %[TRI,DT,LightHandle]=Plot_sbB(CtrlVar,MUA,s,b,B,TRI,DT,AspectRatio,ViewAndLight,LightHandle,sCol,bCol,BCol);
    AspectRatio=1;
    figure ; Plot_sbB(CtrlVar,MUA,F.s,F.b,F.B,[],[],AspectRatio) ; title('F.s, F.b and F.B')
end



if contains(CtrlVar.Inverse.InvertFor,'-B-')


    cbar=UaPlots(CtrlVar,MUA,F,InvFinalValues.B,FigureTitle="B final");
    title('InvFinalValues.B') ;
    subtitle("")
    title(cbar, '(m)');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    colormap(othercolor("Mdarkterrain",32))

    cbar=UaPlots(CtrlVar,MUA,F,InvStartValues.B,FigureTitle="B start");
    title('Bstart')
    subtitle("")
    title(cbar, '(m)')
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    colormap(othercolor("Mdarkterrain",32))

    cbar=UaPlots(CtrlVar,MUA,F,InvFinalValues.B-InvStartValues.B,FigureTitle="B final - B start");
    title('InvFinalValues.B-Bstart') ;
    subtitle("")
    title(cbar, '(m)');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    CM=cmocean('balanced',25,'pivot',0) ; colormap(CM);

    AspectRatio=1;
    figsbB=FindOrCreateFigure("sbB");  clf(figsbB)
    Plot_sbB(CtrlVar,MUA,[],[],F.B,[],[],AspectRatio) ;
    title('B')
    subtitle("")
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);


end






%% Basal drag

[~,~,tb] = CalcBasalTraction(CtrlVar,UserVar,MUA,F) ;
tb(tb<eps)=nan ;
cbar=UaPlots(CtrlVar,MUA,F,tb) ;
title('Basal drag, $\Vert \mathbf{t}_b \Vert$ ','interpreter','latex') ;
title(cbar, '($\mathrm{kPa}$)','interpreter','latex');
subtitle("")
set(gca,'ColorScale','log')
CL=clim ; clim([max(CL(1),0.01) CL(2)])
CM=cmocean('balanced') ; colormap(CM);
%%
% uAdjoint vAdjoint
if isprop(InvFinalValues,'uAdjoint')
    if ~isempty(InvFinalValues.uAdjoint)
        figAV=FindOrCreateFigure("Adjoint variables") ; clf(figAV)
        T=tiledlayout("flow");

        TuAdjoint=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,InvFinalValues.uAdjoint,CreateNewFigure=false);
        title("$u$ Adjoint variable",Interpreter="latex")
        subtitle("")
        CLuAdjoint=clim;

        TvAdjoint=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,InvFinalValues.vAdjoint,CreateNewFigure=false);
        title("$v$ Adjoint variable",Interpreter="latex")
        subtitle("")
        CLvAdjoint=clim;


        axis(TuAdjoint); clim(CLuAdjoint)  ; CM=cmocean('-balanced',25,'pivot',0) ; colormap(TuAdjoint,CM);
        axis(TvAdjoint); clim(CLvAdjoint)  ; CM=cmocean('-balanced',25,'pivot',0) ; colormap(TvAdjoint,CM);
        T.Padding="compact";   T.TileSpacing="tight";

        cbar=UaPlots(CtrlVar,MUA,F,[InvFinalValues.uAdjoint, InvFinalValues.vAdjoint],FigureTitle="Adjoint velocities") ;
        title(cbar,"") ;
        subtitle("") ;
        title("Adjoint velocities")


    end
end
%% Plot velocities and velocity residuals
CtrlVar.VelPlotIntervalSpacing='log10';
if ~exist('x','var')
    x=MUA.coordinates(:,1);
    y=MUA.coordinates(:,2);
end

if ~exist('us','var')
    us=F.ub+F.ud; vs=F.vb+F.vd;
end

if ~exist('usError','var')
    usError=sqrt(spdiags(Meas.usCov));
    vsError=sqrt(spdiags(Meas.vsCov));
    % wsError=sqrt(spdiags(Meas.wsCov));
end
if ~exist('GLgeo','var')
    GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,F.GF,CtrlVar); xGL=[] ; yGL=[] ;
end

%%

figSM=FindOrCreateFigure("Speed misfit") ; clf(figSM)
speedMeas=sqrt(Meas.us.^2+Meas.vs.^2);
speedCalc=sqrt(F.ub.^2+F.vb.^2) ;
ErrSpeed=sqrt(usError.^2+vsError.^2);

T=tiledlayout("flow");

TS1=nexttile;
cbar=UaPlots(CtrlVar,MUA,F,speedMeas,CreateNewFigure=false) ; title('Measured speed') ; set(gca,'ColorScale','log')
title(cbar,"$\|\mathbf{v}_\mathrm{Meas}\|$",interpreter="latex")
subtitle("")
CL=clim;

Ts2=nexttile;
cbar=UaPlots(CtrlVar,MUA,F,speedCalc,CreateNewFigure=false) ; title('Modelled speed') ; set(gca,'ColorScale','log')
title(cbar,"$\|\mathbf{v}_\mathrm{Modelled}\|$",interpreter="latex")
clim(CL);
subtitle("(Same colorbar scale as for measured speed)")

TS3=nexttile;
cbar=UaPlots(CtrlVar,MUA,F,ErrSpeed,CreateNewFigure=false) ; title('Speed measurement error') ; set(gca,'ColorScale','log')
title(cbar,"error",interpreter="latex")
subtitle("")

TS4=nexttile;
D=speedMeas-speedCalc ;
cbar=UaPlots(CtrlVar,MUA,F,D,CreateNewFigure=false) ; title('Measured speed - modelled speed') ; set(gca,'ColorScale','lin')
title(cbar,"$\|\mathbf{v}_\mathrm{Meas}\|-\|\mathbf{v}_{\mathrm{Modelled}}\|$",interpreter="latex")
subtitle("")
T.Padding="tight";   T.TileSpacing="tight";

axis(TS4);
CL=clim;
if CL(1)< 0 && CL(2) > 0
    CM=cmocean('-balanced',25,'pivot',0) ; colormap(TS4,CM);
end

%%
figVmis=FindOrCreateFigure("velocity misfit") ; clf(figVmis)


T=tiledlayout("flow");

nexttile

QuiverColorGHG(x,y,(us-Meas.us)./usError,(vs-Meas.vs)./vsError,CtrlVar);
title("($\mathbf{v}_{\mathrm{modelled}}-\mathbf{v}_{\mathrm{measured}})./\mathbf{v}_{\mathrm{error}}$",Interpreter="latex")
hold on
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
PlotMuaBoundary(CtrlVar,MUA,'b')  ;
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)

nexttile

QuiverColorGHG(x,y,us-Meas.us,vs-Meas.vs,CtrlVar); axis equal ;
title("$\mathbf{v}_{\mathrm{modelled}}-\mathbf{v}_{\mathrm{measured}}$",Interpreter="latex") ;
hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
PlotMuaBoundary(CtrlVar,MUA,'b')  ; xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)



nexttile

[~,~,QuiverPar]=QuiverColorGHG(x,y,Meas.us,Meas.vs,CtrlVar); axis equal ;
title("$\mathbf{v}_{\mathrm{measured}}$",Interpreter="latex")
hold on ;
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
PlotMuaBoundary(CtrlVar,MUA,'b')  ;
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)

nexttile

QuiverPar.QuiverSameVelocityScalingsAsBefore=1;
QuiverColorGHG(x,y,us,vs,QuiverPar); axis equal ;
title("$\mathbf{v}_{\mathrm{modelled}}$",Interpreter="latex")
hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
PlotMuaBoundary(CtrlVar,MUA,'b')  ;
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)
QuiverPar.QuiverSameVelocityScalingsAsBefore=0;

T.Padding="tight";   T.TileSpacing="tight";
%%

if ~isempty(Meas.dhdt)  && contains(CtrlVar.Inverse.Measurements,"-dhdt")

    %%
    figdhdt=FindOrCreateFigure("dh/dt misfit") ; clf(figdhdt)

    T=tiledlayout("flow");

    axdhdt1=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,Meas.dhdt,CreateNewFigure=false);
    title("$\dot{h}_\mathrm{Measured}$",Interpreter="latex") ;
    title(cbar,"$\dot{h}_\mathrm{Measured}$",interpreter="latex")
    subtitle("")
    hold on ;
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)
    CM=cmocean('-balanced',25,'pivot',0) ; colormap(axdhdt1,CM);
    CLmeas=clim;

    axdhdt2=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,F.dhdt,CreateNewFigure=false);
    title("$\dot{h}_{\mathrm{Modelled}}$",Interpreter="latex") ;
    title(cbar,"$\dot{h}_\mathrm{Modelled}$",interpreter="latex")
    subtitle("")
    hold on ;
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)
    CLmod=clim;
    CM=cmocean('-balanced',25,'pivot',0) ; colormap(axdhdt2,CM);

    axdhdt3=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,(F.dhdt-Meas.dhdt)./dhdtError,CreateNewFigure=false);
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)
    title("$\dot{h}_{\mathrm{Modelled}}-\dot{h}_{\mathrm{Measured}}$",Interpreter="latex") ;
    subtitle("")
    CLdiff=clim;
    CM=cmocean('-balanced',25,'pivot',0) ; colormap(axdhdt3,CM);
    title(cbar,"$\Delta \dot{h}$",interpreter="latex")

    axis(axdhdt1); clim(CLmeas)  ; CM=cmocean('-balanced',25,'pivot',0) ; colormap(axdhdt1,CM);
    axis(axdhdt2); clim(CLmod)   ; CM=cmocean('-balanced',25,'pivot',0) ; colormap(axdhdt2,CM);
    axis(axdhdt3); clim(CLdiff)  ; CM=cmocean('-balanced',25,'pivot',0) ; colormap(axdhdt3,CM);
    %%
    T.Padding="tight";   T.TileSpacing="tight";
end



%%
figVel=FindOrCreateFigure("Modelled velocities") ; clf(figVel)
PlotBoundary(MUA.Boundary,MUA.connectivity,MUA.coordinates,CtrlVar,'k')
hold on
QuiverPar.QuiverColorSpeedLimits=[];
QuiverPar.QuiverSameVelocityScalingsAsBefore=0;
QuiverColorGHG(x,y,us,vs,QuiverPar); axis equal ;
title("Modelled horizontal velocities") ;
hold on ;
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,"r");
PlotCalvingFronts(CtrlVar,MUA,F,"b");
%%



UaPlots(CtrlVar,MUA,F,F.dhdt,FigureTitle="dh/dt modelled")
title('Modelled $dh/dt$ (assuming plug flow)','interpreter','latex') ;
subtitle("")
CL=clim;
if CL(1) < 0 && CL(2)>0
    CM=cmocean('balanced',25,'pivot',0) ;
else
    CM=cmocean('balanced',25) ;
end
colormap(CM);
%%  Prior

if contains(CtrlVar.Inverse.InvertFor,"-AGlen-")
    if isscalar(Priors.AGlen)
        Priors.AGlen=Priors.AGlen+zeros(MUA.Nnodes,1);
    end

    cbar=UaPlots(CtrlVar,MUA,F,Priors.AGlen,FigureTitle="APrior") ;
    set(gca,'ColorScale','log')
    title(cbar, '($\mathrm{yr}^{-1}\,\mathrm{kPa}^{-n}$)','interpreter','latex');
    title("$A_{\mathrm{Prior}}$",Interpreter="latex")
    subtitle("")
end

if contains(CtrlVar.Inverse.InvertFor,"-C-")
    cbar=UaPlots(CtrlVar,MUA,F,Priors.C,FigureTitle="CPrior") ;
    set(gca,'ColorScale','log')
    title(cbar, '($\mathrm{m}\,\mathrm{yr}^{-1}\,\mathrm{kPa}^{-m}$)','interpreter','latex');
    title("$C_{\mathrm{Prior}}$",Interpreter="latex")
    subtitle("")
end

if contains(CtrlVar.Inverse.InvertFor,"-B-")
    cbar=UaPlots(CtrlVar,MUA,F,Priors.B,FigureTitle="BPrior") ;
    title(cbar, '($\mathrm{m}$)','interpreter','latex');
    title("$B_{\mathrm{Prior}}$",Interpreter="latex")
    subtitle("")
end



%%

if CtrlVar.Inverse.TestAdjoint.isTrue



    if ~isempty(InvFinalValues.dJdAGlenTest)
        IA=find(~isnan(InvFinalValues.dJdAGlenTest)) ;
        fprintf('------------------------------------ dJ/dA  ---------------------------------------------------------------------\n')
        fprintf('#Node/Ele  dJdA          dJdATest      dJdA-dJdATest     dJdA/dtdATest  (dJdA-dJdATest)/dJdA \n')

        for ii=1:numel(IA)
            I=IA(ii);
            fprintf('%i %15g %15g  %15g  %15g %15g \n',I,...
                InvFinalValues.dJdAGlen(I),...
                InvFinalValues.dJdAGlenTest(I),...
                InvFinalValues.dJdAGlen(I)-InvFinalValues.dJdAGlenTest(I),...
                InvFinalValues.dJdAGlen(I)/InvFinalValues.dJdAGlenTest(I),...
                (InvFinalValues.dJdAGlen(I)-InvFinalValues.dJdAGlenTest(I))/InvFinalValues.dJdAGlen(I))
        end

        figAgrad=FindOrCreateFigure("dJ/dA test") ;  clf(figAgrad)
        plot(InvFinalValues.dJdAGlen,InvFinalValues.dJdAGlenTest,"or") ;
        hold on
        plot(InvFinalValues.dJdAGlen,InvFinalValues.dJdAGlen,"--k") ;

        xlabel("Adjoint $dJ/dA$",Interpreter="latex")  ;
        ylabel("Finite difference $dJ/dA$",Interpreter="latex")
        ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
        axis on ; axis equal tight ; box off
        axis([min(InvFinalValues.dJACTest) max(InvFinalValues.dJACTest) min(InvFinalValues.dJdATest) max(InvFinalValues.dJdATest)])
        title("Comparision betweenadjoint and finite-differences gradient calculations")
        set(gcf,'Color','white')



    end

    if ~isempty(InvFinalValues.dJdCTest)

        IC=find(~isnan(InvFinalValues.dJdCTest)) ;

        fprintf('--------------------------------------- dJ/dC ----------------------------------------------------------------------\n')

        fprintf('#Node/Ele  dJdC          dJdCTest      dJdC-dJdCTest     dJdC/dtdCTest   (dJdC-dJdCTest)/dJdC\n')

        for ii=1:numel(IC)
            I=IC(ii);
            fprintf('%i %15g %15g  %15g  %15g %15g \n',I,...
                InvFinalValues.dJdC(I),...
                InvFinalValues.dJdCTest(I),...
                InvFinalValues.dJdC(I)-InvFinalValues.dJdCTest(I),...
                InvFinalValues.dJdC(I)/InvFinalValues.dJdCTest(I),...
                (InvFinalValues.dJdC(I)-InvFinalValues.dJdCTest(I))/InvFinalValues.dJdC(I))
        end

        %%

        figCgrad=FindOrCreateFigure("dJ/dC test") ;  clf(figCgrad)
        plot(InvFinalValues.dJdC,InvFinalValues.dJdCTest,"or") ;
        hold on
        plot(InvFinalValues.dJdC,InvFinalValues.dJdC,"--k") ;
        xlabel("Adjoint $dJ/dC$",Interpreter="latex")  ;
        ylabel("Finite difference $dJ/dC$",Interpreter="latex")
        ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
        axis on ; axis equal tight;
        axis([min(InvFinalValues.dJdCTest) max(InvFinalValues.dJdCTest) min(InvFinalValues.dJdCTest) max(InvFinalValues.dJdCTest)])
        box off
        title("Comparision betweenadjoint and finite-differences gradient calculations")
        set(gcf,'Color','white')
    end

    if ~isempty(InvFinalValues.dJdBTest)
        fprintf('--------------------------------------- dJ/dB ----------------------------------------------------------------------\n')

        IB=find(~isnan(InvFinalValues.dJdBTest)) ;

        fprintf('#Node/Ele  dJdB          dJdBTest      dJdB-dJdBTest     dJdB/dtdBTest   \n')

        for ii=1:numel(IB)
            I=IB(ii);
            fprintf('%i %15g %15g  %15g  %15g %15g \n',I,...
                InvFinalValues.dJdB(I),...
                InvFinalValues.dJdBTest(I),...
                InvFinalValues.dJdB(I)-InvFinalValues.dJdBTest(I),...
                InvFinalValues.dJdB(I)/InvFinalValues.dJdBTest(I),...
                (InvFinalValues.dJdB(I)-InvFinalValues.dJdBTest(I))/InvFinalValues.dJdB(I))
        end


        fprintf('--------------------------------------------------------------------------------------------------------------------------\n')

        %%

        figBgrad=FindOrCreateFigure("dJ/dB test") ;  clf(figBgrad)
        plot(InvFinalValues.dJdB,InvFinalValues.dJdBTest,"or") ;
        hold on
        plot(InvFinalValues.dJdB,InvFinalValues.dJdB,"--k") ;
        axis equal tight ;
        axis([min(InvFinalValues.dJdBTest) max(InvFinalValues.dJdBTest) min(InvFinalValues.dJdBTest) max(InvFinalValues.dJdBTest)])
        box off
        xlabel("Adjoint $dJ/dB$",Interpreter="latex")  ;
        ylabel("Finite difference $dJ/dB$",Interpreter="latex")
        ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
        title("Comparision between adjoint and finite-differences gradient calculations")
        set(gcf,'Color','white')
    end
    %%
    %[dJdp(iRange) dJdpTest(iRange)   dJdp(iRange)-dJdpTest(iRange) dJdp(iRange)./dJdpTest(iRange)]
    % iRange=find(~isnan(dJdpTest));
    % dJdpTest(iRange)-dJdp(iRange))/norm(dJdp(iRange)))
    % fprintf('Norm test: ||dJdpTest-dJdp||/||dJdp||= %g \n ',norm(dJdpTest(iRange)-dJdp(iRange))/norm(dJdp(iRange)))

    %%

    if ~(isempty(InvFinalValues.dJdC) && isempty(InvFinalValues.dJdCTest))

        IFigC=FindOrCreateFigure("dJ/dC test over mesh") ; clf(IFigC);

        TileC=tiledlayout("flow");

        TdJdC=nexttile;

        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdC) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('dJdC Adjoint directional derivative')

        TdJdCTest=nexttile;
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdCTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('$dJ/dC$ Brute force directional derivative','interpreter','latex')


        TdJdCDiff=nexttile;
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdC-InvFinalValues.dJdCTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('Difference between adjoint and brute force directional derivatives')


        TdJdCRatio=nexttile;
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdC./InvFinalValues.dJdCTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('Ratio between adjoint and brute force directional derivatives')

        TileC.TileSpacing='tight';
        TileC.Padding='tight';
        %IFigC.Position=[948.43 41.571 1246.3 1115.4];
        %%
    end


    if ~(isempty(InvFinalValues.dJdAGlen) && isempty(InvFinalValues.dJdAGlenTest))


        IFigAGlen=FindOrCreateFigure("dJ/dA test over mesh") ; clf(IFigAGlen);

        TileA=tiledlayout("flow");

        nexttile
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlen) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('dJdAGlen Adjoint directional derivative')

        nexttile
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlenTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('dJdAGlen Brute force directional derivative')


        nexttile
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlen-InvFinalValues.dJdAGlenTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('Difference between adjoint and brute force directional derivatives')

        nexttile
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlen./InvFinalValues.dJdAGlenTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('Ratio between adjoint and brute force directional derivatives')

        TileA.TileSpacing='tight';
        TileA.Padding='tight';
        % IFigAGlen.Position=[1.5714 41.571 1096 1115.4];
        %%
    end


    if ~(isempty(InvFinalValues.dJdB) && isempty(InvFinalValues.dJdBTest))

        %%

        IFigb=FindOrCreateFigure("dJ/dB test over mesh") ; clf(IFigb)
        TileB=tiledlayout("flow");
        TdJdB=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,InvFinalValues.dJdB,PlotUnderMesh=true,CreateNewFigure=false);
        title("$dJ/dB$ Adjoint directional derivative")
        subtitle("")

        TdJdBTest=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,InvFinalValues.dJdBTest,PlotUnderMesh=true,CreateNewFigure=false);
        title('$dJ/dB$ Brute force directional derivative',Interpreter='latex')

        TdJdBDiff=nexttile;
        cbar=UaPlots(CtrlVar,MUA,F,InvFinalValues.dJdB-InvFinalValues.dJdBTest,PlotUnderMesh=true,CreateNewFigure=false);
        title('Difference between adjoint and brute force derivatives')
        subtitle("")
        axis(TdJdBDiff) ; CM=cmocean('balanced',25,'pivot',0) ; colormap(TdJdBDiff,CM);

        TdJdBRatio=nexttile;
        UaPlots(CtrlVar,MUA,F,InvFinalValues.dJdB./InvFinalValues.dJdBTest,PlotUnderMesh=true,CreateNewFigure=false) ;
        title('Ratio between adjoint and brute force derivatives')
        subtitle("")

        %  IFigb.Position=[1.5714 41.571 1096 1115.4];
        TileB.TileSpacing='tight';
        TileB.Padding='tight';
        %CL=TdJdB.CLim; TdJdBTest.CLim=CL;  TdJdBDiff.CLim=CL; TdJdBRatio.CLim=[0.7 1.3];
        axis(TdJdBDiff) ; CM=cmocean('balanced',25,'pivot',0) ; colormap(TdJdBDiff,CM);
        %%
    end

else

    if ~isempty(InvFinalValues.dJdAGlen)

        PM=CtrlVar.Inverse.AdjointGradientPreMultiplier;
        figPMdJdA=FindOrCreateFigure(PM+"\dJdA "); clf(figPMdJdA) ;
        UaPlots(CtrlVar,MUA,F,InvFinalValues.dJdAGlen,CreateNewFigure=false);

        if PM=="M"
            T="$\nabla_A J = M^{-1} dJ/dA$";
        else
            T="$\nabla_A J=dJ/dA$";
        end
        title(T,Interpreter="latex")

        subtitle("")
        cl=clim;
        if min(cl) <0 && max(cl)> 0
            CM=cmocean('balanced',25,'pivot',0) ; colormap(figPMdJdA,CM);
        else
            CM=cmocean('balanced',25) ; colormap(figPMdJdA,CM);
        end

        if PM=="I"

            if isempty(MUA.M)
                MUA.M=MassMatrix2D1dof(MUA);
            end


            dJdA=MUA.M\InvFinalValues.dJdAGlen;
            figMdJdA=FindOrCreateFigure("M\dJdA "+PM); clf(figMdJdA)
            UaPlots(CtrlVar,MUA,F,dJdA,CreateNewFigure=false);
            T="$\nabla_C J =M^{-1} dJdA$" ;
            title(T,interpreter="latex")
            subtitle("")
            cl=clim;
            if min(cl) <0 && max(cl)> 0
                CM=cmocean('balanced',25,'pivot',0) ; colormap(figMdJdA,CM);
            else
                CM=cmocean('balanced',25) ; colormap(figMdJdA,CM);
            end
        end

    end



    if ~isempty(InvFinalValues.dJdC)

        PM=CtrlVar.Inverse.AdjointGradientPreMultiplier;
        figPMdJdC=FindOrCreateFigure(PM+"\dJdC "); clf(figPMdJdC)
        UaPlots(CtrlVar,MUA,F,InvFinalValues.dJdC,CreateNewFigure=false);
        if PM=="M"
            T="$\nabla_C J = M^{-1} dJ/dC$";
        else
            T="$\nabla_C J=dJ/dC$";
        end

        title(T,Interpreter="latex")
        subtitle("")
        cl=clim;
        if min(cl) <0 && max(cl)> 0
            CM=cmocean('balanced',25,'pivot',0) ; colormap(figPMdJdC,CM);
        else
            CM=cmocean('balanced',25) ; colormap(figPMdJdC,CM);
        end

        if PM=="I"
            dJdC=MUA.M\InvFinalValues.dJdC;
            figMdJdC=FindOrCreateFigure("M\dJdC "+PM); clf(figMdJdC)
            UaPlots(CtrlVar,MUA,F,dJdC,CreateNewFigure=false);
            T="$\nabla_C J =M^{-1} dJ/dC$";
            title(T,interpreter="latex")
            subtitle("")
            cl=clim;
            if min(cl) <0 && max(cl)> 0
                CM=cmocean('balanced',25,'pivot',0) ; colormap(figMdJdC,CM);
            else
                CM=cmocean('balanced',25) ; colormap(figMdJdC,CM);
            end
        end
    end

    if ~isempty(InvFinalValues.dJdB)

        PM=CtrlVar.Inverse.AdjointGradientPreMultiplier;
        figdJdB=FindOrCreateFigure("dJdB"+PM); clf(figdJdB)
        UaPlots(CtrlVar,MUA,F,InvFinalValues.dJdB,CreateNewFigure=false);
        if PM=="M"
            T="$\nabla_B J = M^{-1} dJ/dB$";
        else
            T="$\nabla_B J=dJ/dB$";
        end
        title(T,Interpreter="latex")
        subtitle("")

        cl=clim;
        if min(cl) <0 && max(cl)> 0
            CM=cmocean('balanced',25,'pivot',0) ; colormap(figdJdB,CM);
        else
            CM=cmocean('balanced',25) ; colormap(figdJdB,CM);
        end

        if PM=="I"
            dJdB=MUA.M\InvFinalValues.dJdB;
            figMB=FindOrCreateFigure("M\dJdB "+PM); clf(figMB)
            UaPlots(CtrlVar,MUA,F,dJdB,CreateNewFigure=false);
            T="$\nabla_B J =M^{-1} dJ/dB$";
            title(T,interpreter="latex")
            subtitle("")
            cl=clim;
            if min(cl) <0 && max(cl)> 0
                CM=cmocean('balanced',25,'pivot',0) ; colormap(figMB,CM);
            else
                CM=cmocean('balanced',25) ; colormap(figMB,CM);
            end
        end


    end


    %subplot(3,1,3) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dRdp) ; title('dRdp')


    %IFigGradients.Position=[1098.7 638.71 1096 518.29];
    %%
    CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;

    if contains(lower(CtrlVar.Inverse.InvertFor),'c') && contains(lower(CtrlVar.Inverse.InvertFor),'aglen')

        figAC=FindOrCreateFigure("True and estimated A and C"); clf(figAC);

        T=tiledlayout("flow");

        TTrueC=nexttile;
        UaPlots(CtrlVar,MUA,F,Priors.TrueC,CreateNewFigure=false) ;
        title("True $C$",interpreter="latex") ; set(gca,'ColorScale','log')
        subtitle("")

        TRetrievedC=nexttile;
        UaPlots(CtrlVar,MUA,F,InvFinalValues.C,CreateNewFigure=false) ;
        title("Retrieved $C$",interpreter="latex") ;
        set(gca,'ColorScale','log')
        subtitle("")


        CCbarLink=linkprop([TTrueC TRetrievedC],'CLim') ; assignin('base','CbarLink_clim',CCbarLink)

        TTrueA=nexttile;
        UaPlots(CtrlVar,MUA,F,Priors.TrueAGlen,CreateNewFigure=false) ;
        title("True $A$",Interpreter="latex") ;
        set(gca,'ColorScale','log')
        subtitle("")

        TRetrievedA=nexttile;
        UaPlots(CtrlVar,MUA,F,InvFinalValues.AGlen,CreateNewFigure=false) ;
        title("Retrieved $A$",Interpreter="latex") ; set(gca,'ColorScale','log')
        subtitle("")


        ACbarLink=linkprop([TTrueA TRetrievedA],'CLim') ; assignin('base','CbarLink_clim',ACbarLink)

        %figC.Position=[400 200 1300 800];
        T.Padding="tight";
        T.TileSpacing="tight";

    end





    if contains(lower(CtrlVar.Inverse.InvertFor),'c')
        if ~isempty(Priors.TrueC)  && ~anynan(Priors.TrueC)


            figC=FindOrCreateFigure("True and estimated C"); clf(figC);

            T=tiledlayout("flow");

            nexttile
            UaPlots(CtrlVar,MUA,F,Priors.TrueC,CreateNewFigure=false) ;
            title("True $C$",interpreter="latex") ; set(gca,'ColorScale','log')
            subtitle("")

            nexttile
            UaPlots(CtrlVar,MUA,F,InvFinalValues.C,CreateNewFigure=false) ;
            title("Retrieved $C$",interpreter="latex") ; set(gca,'ColorScale','log')
            subtitle("")

            nexttile

            D=abs(Priors.TrueC-InvFinalValues.C) ;
            cbar=UaPlots(CtrlVar,MUA,F,D,CreateNewFigure=false) ;
            title('abs(True C - Retrieved C)') ; set(gca,'ColorScale','log')
            title(cbar,"$|C-\tilde{C}|$",interpreter="latex")
            subtitle("")

            nexttile
            UaPlots(CtrlVar,MUA,F,InvStartValues.C,CreateNewFigure=false);
            title("C at start of inversion") ; set(gca,'ColorScale','log')
            subtitle("")

            nexttile
            UaPlots(CtrlVar,MUA,F,Priors.C,CreateNewFigure=false) ;
            title('Prior C') ; set(gca,'ColorScale','log')
            subtitle("")


            nexttile
            UaPlots(CtrlVar,MUA,F,InvFinalValues.C-Priors.C,CreateNewFigure=false);
            title("Retrieved C -  Prior C ") ; set(gca,'ColorScale','log')
            subtitle("")

            %figC.Position=[400 200 1300 800];
            T.Padding="tight";
            T.TileSpacing="tight";


        end
    end

    if contains(lower(CtrlVar.Inverse.InvertFor),'aglen')

        if ~isempty(Priors.TrueAGlen) && ~anynan(Priors.TrueAGlen)



            tFig1=FindOrCreateFigure("True and estimated AGlen"); clf(tFig1 );

            T=tiledlayout("flow");

            nexttile
            UaPlots(CtrlVar,MUA,F,Priors.TrueAGlen,CreateNewFigure=false) ;
            title('True AGlen') ; set(gca,'ColorScale','log')
            subtitle("")

            nexttile
            UaPlots(CtrlVar,MUA,F,InvFinalValues.AGlen,CreateNewFigure=false) ;
            title('Retrieved AGlen') ; set(gca,'ColorScale','log')
            subtitle("")

            nexttile

            D=abs(Priors.TrueAGlen-InvFinalValues.AGlen) ;
            cbar=UaPlots(CtrlVar,MUA,F,D,CreateNewFigure=false) ;
            title('abs(True A -Retrieved A)') ; set(gca,'ColorScale','log')
            title(cbar,"$|A-\tilde{A}|$",interpreter="latex")
            subtitle("")

            nexttile
            UaPlots(CtrlVar,MUA,F,Priors.AGlen,CreateNewFigure=false) ;
            title('Prior AGlen') ; set(gca,'ColorScale','log')
            subtitle("")

            T.Padding="tight";   T.TileSpacing="tight";



        end
    end


    if contains(CtrlVar.Inverse.InvertFor,'-B-')

        %% B

        PlotBedrockInversionFields(CtrlVar,MUA,F,Priors,InvFinalValues,InvStartValues,Meas)


    end



    %%
    if ~isempty(RunInfo.Inverse.J)

        figIPO=FindOrCreateFigure("Inverse Parameter Optimisation");
        clf(figIPO)
        hold off
        yyaxis left
        semilogy(RunInfo.Inverse.Iterations,RunInfo.Inverse.J,'-bo','LineWidth',2)
        ylabel("$J$",'interpreter','latex')

        if ~isempty(RunInfo.Inverse.GradNorm)  && ~all(isnan(RunInfo.Inverse.GradNorm)) ...
                &&  numel(RunInfo.Inverse.Iterations) == numel(RunInfo.Inverse.GradNorm)

            hold off
            yyaxis right
            semilogy(RunInfo.Inverse.Iterations,RunInfo.Inverse.GradNorm,'-r+')
            ylabel('Norm of gradient ($l_2$)','interpreter','latex')
            legend('Objective function','$\| \hbox{Gradient} \|$','Location','northeast','interpreter','latex')

        end

        yyaxis left
        xlabel('Inverse iteration','interpreter','latex');
        hold off

        if ~all(isnan(RunInfo.Inverse.R))

            figJIR=FindOrCreateFigure('J=I+R');
            clf(figJIR)
            hold off
            yyaxis left
            semilogy(RunInfo.Inverse.Iterations,RunInfo.Inverse.J,'-bo','LineWidth',2)
            ylabel('$J$','interpreter','latex')

            hold on
            semilogy(RunInfo.Inverse.Iterations,RunInfo.Inverse.I,'-gx')
            ylabel('$J$ and $I$',Interpreter='latex')

            yyaxis right
            semilogy(RunInfo.Inverse.Iterations,RunInfo.Inverse.R,'-r+')
            ylabel('$R$','interpreter','latex')
            xlabel('Inverse iteration','interpreter','latex');
            legend('Objective function','$I$','$R$','Location','southwest','interpreter','latex')

        end

        yyaxis left
        xlabel('Inverse iteration','interpreter','latex') ;
        hold off
    end


end

FindOrCreateFigure("Forward Boundary Conditions") ;
PlotBoundaryConditions(CtrlVar,MUA,BCs)

fprintf("...done \n")

FindOrCreateFigure("Adjoint Boundary Conditions") ;
CtrlVar.BCsType="adjoint";

PlotBoundaryConditions(CtrlVar,MUA,BCsAdjoint);

fprintf("...done \n")

end
