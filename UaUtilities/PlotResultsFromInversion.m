




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

    ColorbarLimits=10.^[mean(log10(InvFinalValues.AGlen))-4*std(log10(InvFinalValues.AGlen))  mean(log10(InvFinalValues.AGlen))+4*std(log10(InvFinalValues.AGlen))];
    if ColorbarLimits(1)==ColorbarLimits(2)
        Eps=10*eps(ColorbarLimits(1));
        ColorbarLimits(1)=ColorbarLimits(1)-Eps;
        ColorbarLimits(2)=ColorbarLimits(2)+Eps;
    end

    figCeI=FindOrCreateFigure("Change in A during inversion") ; clf(figCeI)


    T=tiledlayout("flow");


    T1=nexttile ;
    UaPlots(CtrlVar,MUA,F,F.AGlen,CreateNewFigure=false,logColorbar=true);
    xlabel(CtrlVar.PlotsXaxisLabel,Interpreter="latex")  ; ylabel(CtrlVar.PlotsYaxisLabel,Interpreter="latex")
    title("$A$ at end of inversion",Interpreter="latex"); subtitle("")
    cbar=colorbar;
    title(cbar, '($\mathrm{m}\,\mathrm{yr}^{-1}\,\mathrm{kPa}^{-m}$)','interpreter','latex');
    CM=cmocean('-ice',15) ; colormap(CM);
    clim(ColorbarLimits)

    T2=nexttile ;
    UaPlots(CtrlVar,MUA,F,InvStartValues.AGlen,CreateNewFigure=false,logColorbar=true);
    xlabel(CtrlVar.PlotsXaxisLabel,Interpreter="latex")  ; ylabel(CtrlVar.PlotsYaxisLabel,Interpreter="latex")
    title("$A$ at start of current inversion run",Interpreter="latex") ; subtitle("")
    cbar=colorbar;
    title(cbar, '($\mathrm{m}\,\mathrm{yr}^{-1}\,\mathrm{kPa}^{-m}$)','interpreter','latex');
    CM=cmocean('-ice',15) ; colormap(CM);
    clim(ColorbarLimits)


    T3=nexttile ;
    dC=log10(InvFinalValues.AGlen)-log10(InvStartValues.AGlen);
    cbar=UaPlots(CtrlVar,MUA,F,dC,CreateNewFigure=false);
    xlabel(CtrlVar.PlotsXaxisLabel,Interpreter="latex")  ; ylabel(CtrlVar.PlotsYaxisLabel,Interpreter="latex")
    title("Change in $A$ during current inversion run",Interpreter="latex") ;
    subtitle("$\log(A_{\mathrm{End}})-\log(A_{\mathrm{Start}})$",Interpreter="latex")
    title(cbar, '($\mathrm{m}\,\mathrm{yr}^{-1}\,\mathrm{kPa}^{-m}$)','interpreter','latex');
    set(figCeI,CurrentAxes=T3) ;
    CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);


    T.Padding="tight";   T.TileSpacing="tight";

end

%%
if contains(upper(CtrlVar.Inverse.InvertFor),'C')

    ColorbarLimits=10.^[mean(log10(InvFinalValues.C))-4*std(log10(InvFinalValues.C))  mean(log10(InvFinalValues.C))+4*std(log10(InvFinalValues.C))];
    if ColorbarLimits(1)==ColorbarLimits(2)
        Eps=10*eps(ColorbarLimits(1));
        ColorbarLimits(1)=ColorbarLimits(1)-Eps;
        ColorbarLimits(2)=ColorbarLimits(2)+Eps;
    end

    figCeI=FindOrCreateFigure("Change in C during inversion") ; clf(figCeI)


    T=tiledlayout("flow");


    T1=nexttile ;
    UaPlots(CtrlVar,MUA,F,F.C,CreateNewFigure=false,logColorbar=true);
    xlabel(CtrlVar.PlotsXaxisLabel,Interpreter="latex")  ; ylabel(CtrlVar.PlotsYaxisLabel,Interpreter="latex")
    title("$C$ at end of inversion",Interpreter="latex"); subtitle("")
    cbar=colorbar;
    title(cbar, '($\mathrm{m}\,\mathrm{yr}^{-1}\,\mathrm{kPa}^{-m}$)','interpreter','latex');
    CM=cmocean('-ice',15) ; colormap(CM);
    clim(ColorbarLimits)

    T2=nexttile ;
    UaPlots(CtrlVar,MUA,F,InvStartValues.C,CreateNewFigure=false,logColorbar=true);
    xlabel(CtrlVar.PlotsXaxisLabel,Interpreter="latex")  ; ylabel(CtrlVar.PlotsYaxisLabel,Interpreter="latex")
    title("$C$ at start of current inversion run",Interpreter="latex") ; subtitle("")
    cbar=colorbar;
    title(cbar, '($\mathrm{m}\,\mathrm{yr}^{-1}\,\mathrm{kPa}^{-m}$)','interpreter','latex');
    CM=cmocean('-ice',15) ; colormap(CM);
    clim(ColorbarLimits)


    T3=nexttile ;
    dC=log10(InvFinalValues.C)-log10(InvStartValues.C);
    cbar=UaPlots(CtrlVar,MUA,F,dC,CreateNewFigure=false);
    xlabel(CtrlVar.PlotsXaxisLabel,Interpreter="latex")  ; ylabel(CtrlVar.PlotsYaxisLabel,Interpreter="latex")
    title("Change in $C$ during current inversion run",Interpreter="latex") ;
    subtitle("$\log(C_{\mathrm{End}})-\log(C_{\mathrm{Start}})$",Interpreter="latex")
    title(cbar, '($\mathrm{m}\,\mathrm{yr}^{-1}\,\mathrm{kPa}^{-m}$)','interpreter','latex');

    set(figCeI,CurrentAxes=T3) ;
    CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);

    T.Padding="tight";   T.TileSpacing="tight";

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



if contains(CtrlVar.Inverse.InvertFor,"-B-")

    figB=FindOrCreateFigure("Change in B during inversion") ; clf(figB)


    T=tiledlayout("flow");


    T1=nexttile ;
    cbar=UaPlots(CtrlVar,MUA,F,InvFinalValues.B,CreateNewFigure=false);
    title("$B$ at end of inversion run",Interpreter="latex")
    subtitle("")
    title(cbar, '(m)');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    colormap(T1,othercolor("Mdarkterrain",32))

    T2=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,InvStartValues.B,CreateNewFigure=false);
    title("$B$ at start of current inversion run",Interpreter="latex")
    subtitle("")
    title(cbar, '(m)')
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    colormap(T2,othercolor("Mdarkterrain",32))

    T3=nexttile ;
    cbar=UaPlots(CtrlVar,MUA,F,InvFinalValues.B-InvStartValues.B,CreateNewFigure=false);
    title("Change in $B$ during current inversion run",Interpreter="latex")
    subtitle("")
    title(cbar, '(m)');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);

    set(figB,CurrentAxes=T3) ; CM=cmocean('balanced',25,'pivot',0) ; colormap(T3,CM);
    set(figB,CurrentAxes=T1) ;  colormap(T1,othercolor("Mdarkterrain",32))
    set(figB,CurrentAxes=T2) ;  colormap(T2,othercolor("Mdarkterrain",32))

    T.Padding="tight";   T.TileSpacing="tight";



    AspectRatio=1;
    figsbB=FindOrCreateFigure("sbB");  clf(figsbB)
    Plot_sbB(CtrlVar,MUA,[],[],F.B,[],[],AspectRatio) ;
    title('B')
    subtitle("")
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);


    if ~isempty(Meas.B)

        figBmeas=FindOrCreateFigure("Direct measurements of B") ; clf(figBmeas)

        T=tiledlayout("flow");

        T1=nexttile ;
        cbar=UaPlots(CtrlVar,MUA,F,Meas.B,CreateNewFigure=false);
        title("Direct $B$ Measurements",Interpreter="latex")
        subtitle("")
        title(cbar, '(m)')
        xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
        colormap(othercolor("Mdarkterrain",32))

        T2=nexttile ;

        Berr=full(sqrt(diag(Meas.BCov)))  ; Berr(~isfinite(Berr))=nan ;
        cbar=UaPlots(CtrlVar,MUA,F,Berr,CreateNewFigure=false,logColorbar=true);

        title("Direct $B$ Measurement Errors",Interpreter="latex")
        subtitle("")
        title(cbar, '(m)')
        xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);

        T.Padding="tight";   T.TileSpacing="tight";

        set(figBmeas,CurrentAxes=T1) ;  colormap(T1,othercolor("Mdarkterrain",32))
        set(figBmeas,CurrentAxes=T2) ;  colormap(T2,parula)


    end


end



%% Basal drag

[~,~,tb] = CalcBasalTraction(CtrlVar,UserVar,MUA,F) ;
tb(tb<eps)=nan ;
cbar=UaPlots(CtrlVar,MUA,F,tb,FigureTitle="Basal drag") ;
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
        axis(TvAdjoint); clim(CLuAdjoint)  ; CM=cmocean('-balanced',25,'pivot',0) ; colormap(TvAdjoint,CM);
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

figSM=FindOrCreateFigure("Misfit:Speed") ; clf(figSM)
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


axis(TS4);
CL=clim;
if CL(1)< 0 && CL(2) > 0
    CM=cmocean('-balanced',25,'pivot',0) ; colormap(TS4,CM);
end

T.Padding="tight";   T.TileSpacing="tight";

%%
figVmis=FindOrCreateFigure("Misfit:Velocities") ; clf(figVmis)


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
    figdhdt=FindOrCreateFigure("Misfit:dh/dt") ; clf(figdhdt)

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
    title("$(\dot{h}_{\mathrm{Modelled}}-\dot{h}_{\mathrm{Measured}})/\dot{h}_{\mathrm{Error}}$",Interpreter="latex") ;
    subtitle("")
    CLdiff=clim;
    CM=cmocean('-balanced',25,'pivot',0) ; colormap(axdhdt3,CM);
    %title(cbar,"$\Delta \dot{h}/\dot{h}_{\mathrm{Error}}$",interpreter="latex")

    axdhdt4=nexttile;

    cbar=UaPlots(CtrlVar,MUA,F,dhdtError,CreateNewFigure=false,logColorbar=true);
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)
    title("$\dot{h}_{\mathrm{Error}}$",Interpreter="latex") ;
    subtitle("")
    title(cbar,"$\dot{h}_{\mathrm{Error}}$",interpreter="latex")

    axis(axdhdt1); clim(CLmeas)  ; CM=cmocean('-balanced',25,'pivot',0) ; colormap(axdhdt1,CM);
    axis(axdhdt2); clim(CLmod)   ; CM=cmocean('-balanced',25,'pivot',0) ; colormap(axdhdt2,CM);
    axis(axdhdt3); clim(CLdiff)  ; CM=cmocean('-balanced',25,'pivot',0) ; colormap(axdhdt3,CM);


    T.Padding="tight";   T.TileSpacing="tight";

    %%

    figflux=FindOrCreateFigure("flux divergence etc") ; clf(figflux)

    T=tiledlayout("flow");

    flux1=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,F.rho.*F.ab,CreateNewFigure=false);
    title(cbar,"$(\mathrm{kg}\,\mathrm{m^{-2}}\,\mathrm{a^{-1}})$",interpreter="latex")
    title("$\rho\,a_b$",interpreter="latex")
    subtitle("")
    hold on ;
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)
    T1clim=clim;
    CM=cmocean('-balanced',25,'pivot',0) ; colormap(flux1,CM);
   


    flux2=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,F.rho.*F.as,CreateNewFigure=false);
    title("$\rho \,\dot{a}_s$",Interpreter="latex") ;
    title(cbar,"$(\mathrm{kg}\,\mathrm{m^{-2}}\,\mathrm{a^{-1}})$",interpreter="latex")
    subtitle("")
    hold on ;
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)
    T2clim=clim;
    CM=cmocean('-balanced',25,'pivot',0) ; colormap(flux2,CM);
    CLmeas=clim;

    qx=F.rho.*F.ub.*F.h ; qy=F.rho.*F.vb.*F.h;
    [dqxdx,dqxdy]=calcFEderivativesMUA(qx,MUA,CtrlVar);
    [dqydx,dqydy]=calcFEderivativesMUA(qy,MUA,CtrlVar);
    qdiv=dqxdx+dqydy;
    qdiv=ProjectFintOntoNodes(CtrlVar,MUA,qdiv) ;

    flux3=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,qdiv,CreateNewFigure=false);
    title("$\nabla q  \; (\mathrm{kg}\,\mathrm{m^{-2}}\,\mathrm{a^{-1}})$",Interpreter="latex") ;
    title(cbar,"")
    subtitle("")
    hold on ;
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)
    T3clim=clim;
    CM=cmocean('-balanced',25,'pivot',0) ; colormap(flux3,CM);

    [dudx,dudy]=calcFEderivativesMUA(F.ub,MUA,CtrlVar);
    [dvdx,dvdy]=calcFEderivativesMUA(F.vb,MUA,CtrlVar);
    vdiv=dudx+dvdy;
    vdiv=ProjectFintOntoNodes(CtrlVar,MUA,vdiv) ;

    flux4=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,vdiv,CreateNewFigure=false);
    title("$\nabla \cdot v  \; (\mathrm{a^{-1}})$",Interpreter="latex") ;
    title(cbar,"")
    subtitle("")
    hold on ;
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)
    T4clim=clim;
    CM=cmocean('-balanced',25,'pivot',0) ; colormap(flux4,CM);

    dhdtEst=(F.rho.*(F.as+F.ab)-qdiv)./F.rho ;

    flux5=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,dhdtEst,CreateNewFigure=false);
    title("$a - (\nabla \cdot q )/\rho  \; (\mathrm{m}\;\mathrm{a^{-1}})$",Interpreter="latex") ;
    title(cbar,"")
    subtitle("")
    hold on ;
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)
    T5clim=clim;
    CM=cmocean('-balanced',25,'pivot',0) ; colormap(flux5,CM);


    [dbdx,dbdy]=calcFEderivativesMUA(F.b,MUA,CtrlVar);
    [dbdx,dbdy]=ProjectFintOntoNodes(CtrlVar,MUA,dbdx,dbdy) ;
    db=sqrt(dbdx.*dbdx+dbdy.*dbdy); 
    flux6=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,db,CreateNewFigure=false,logColorbar=true);
    title("Norm of lower ice surface gradients $\| \nabla b \|$",Interpreter="latex") ;
    title(cbar,"")
    subtitle("")
    hold on ;
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)
 
  
    [dsdx,dsdy]=calcFEderivativesMUA(F.s,MUA,CtrlVar);
    [dsdx,dsdy]=ProjectFintOntoNodes(CtrlVar,MUA,dsdx,dsdy) ;
    ds=sqrt(dsdx.*dsdx+dsdy.*dsdy); 
    flux7=nexttile;
    cbar=UaPlots(CtrlVar,MUA,F,ds,CreateNewFigure=false,logColorbar=true);
    title("Norm of upper ice surface gradients $\| \nabla s \|$",Interpreter="latex") ;
    title(cbar,"")
    subtitle("")
    hold on ;
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)



    axis(flux1); clim(T1clim) ; CM=cmocean('-balanced',25,'pivot',0) ; colormap(flux1,CM);
    axis(flux2); clim(T2clim) ; CM=cmocean('-ice',25) ; colormap(flux2,CM);
    axis(flux3); clim(T3clim) ; CM=cmocean('-balanced',25,'pivot',0) ; colormap(flux3,CM);
    axis(flux4); clim(T4clim) ; CM=cmocean('-balanced',25,'pivot',0) ; colormap(flux4,CM);
    axis(flux5); clim(T5clim) ; CM=cmocean('-balanced',25,'pivot',0) ; colormap(flux5,CM);

    T.Padding="tight";   T.TileSpacing="tight";

    %%


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

    PlotTestAdjoint(CtrlVar,MUA,F,InvFinalValues)


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

            if isempty(MUA.M)
                MUA.M=MassMatrix2D1dof(MUA);
            end


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
        figdJdB=FindOrCreateFigure(PM+"\dJdB"); clf(figdJdB)
        UaPlots(CtrlVar,MUA,F,InvFinalValues.dJdB,CreateNewFigure=false);
        if PM=="M"
            T="($\nabla_B J, \delta B)_{L^2} =  \delta_B J[\delta B]$";
        elseif PM=="H1"
            T="($\nabla_B J, \delta B)_{H^1} =  \delta_B J[\delta B]$";
        else
            T="($\nabla_B J, \delta B)_{l^2} =  \delta_B J[\delta B]$";
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
            if isempty(MUA.M)
                MUA.M=MassMatrix2D1dof(MUA);
            end
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
