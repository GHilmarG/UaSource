function PlotResultsFromInversion(UserVar,CtrlVar,MUA,BCs,F,~,~,InvStartValues,InvFinalValues,Priors,Meas,~,RunInfo)


%%
%
% PlotResultsFromInversion(UserVar,CtrlVar,MUA,BCs,F,~,~,InvStartValues,InvFinalValues,Priors,Meas,~,RunInfo)
%
% Does what it says on the tin.
%
%  Example:
%
% load InversionRestartFile
% PlotResultsFromInversion(UserVar,CtrlVar,MUA,BCs,F,~,~,InvStartValues,InvFinalValues,Priors,Meas,~,RunInfo)
%
% It is also possible to enter the name of the restart file as the first, and only, argument. Then the restart file will be
% first loaded, and then plotted.
%
%
% Note: This function is used by Ua for plotting results from an inversion.
%
%%

if isstring(UserVar) && isfile(UserVar)

    fprintf("loading and plotting results from %s \n",UserVar)

    load(UserVar,"UserVarInRestartFile","CtrlVarInRestartFile","MUA","BCs","F","InvStartValues","InvFinalValues","Priors","Meas","RunInfo") ;

    CtrlVar=CtrlVarInRestartFile;
    UserVar=UserVarInRestartFile;

end

%%

fprintf(' Plotting results from inversion...')

CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
us=F.ub+F.ud; vs=F.vb+F.vd;

x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);

GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,F.GF,CtrlVar); xGL=[] ; yGL=[] ;

%%

if ~isempty(Meas.dhdt)
    Iplot=2 ; Jplot=3;
else
    Iplot=2 ; Jplot=2;
end
Kplot=0;


fig=FindOrCreateFigure('Measurements') ; clf(fig)

T=tiledlayout("flow");

nexttile

cbar=UaPlots(CtrlVar,MUA,F,Meas.us,CreateNewFigure=false);

xlabel(CtrlVar.PlotsXaxisLabel,'interpreter','latex');
ylabel(CtrlVar.PlotsYaxisLabel,'interpreter','latex');
title('us Meas on numerical grid') ;


nexttile

cbar=UaPlots(CtrlVar,MUA,F,Meas.vs,CreateNewFigure=false);

xlabel(CtrlVar.PlotsXaxisLabel,'interpreter','latex');
ylabel(CtrlVar.PlotsYaxisLabel,'interpreter','latex');
title('vs Meas on numerical grid') ;

if ~isempty(Meas.dhdt)  && contains(CtrlVar.Inverse.Measurements,"-dhdt")
    
    nexttile
    cbar=UaPlots(CtrlVar,MUA,F,Meas.dhdt,CreateNewFigure=false);

    xlabel(CtrlVar.PlotsXaxisLabel,'interpreter','latex');
    ylabel(CtrlVar.PlotsYaxisLabel,'interpreter','latex');
    title('dh/dt Meas on numerical grid') ;
end

usError=sqrt(spdiags(Meas.usCov));
vsError=sqrt(spdiags(Meas.vsCov));
dhdtError=sqrt(spdiags(Meas.dhdtCov));



nexttile
cbar=UaPlots(CtrlVar,MUA,F,usError,CreateNewFigure=false);
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
title('us error on numerical grid') ;


nexttile
cbar=UaPlots(CtrlVar,MUA,F,vsError,CreateNewFigure=false);
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
title('vs error on numerical grid') ;

if ~isempty(Meas.dhdt)  && contains(CtrlVar.Inverse.Measurements,"-dhdt")
  
    nexttile
    cbar=UaPlots(CtrlVar,MUA,F,dhdtError,CreateNewFigure=false);
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    title('dh/dt error on numerical grid') ;

end
T.Padding="tight";   T.TileSpacing="tight";


%%
if contains(upper(CtrlVar.Inverse.InvertFor),'A')

    fig=FindOrCreateFigure('A at the end of inversion') ; clf(fig)
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

    fig=FindOrCreateFigure('A at the start of inversion') ; clf(fig)
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

    fig=FindOrCreateFigure('Change in A during inversion run') ; clf(fig)
    PlotMeshScalarVariable(CtrlVar,MUA,log10(InvFinalValues.AGlen)-log10(InvStartValues.AGlen));
    title('log10(InvFinalValues.AGlen)-log10(InvStartValues.AGlen)') ; 
    cbar=colorbar; title(cbar, '($\mathrm{a}^{-1}$ $\mathrm{kPa}^{-3}$)',interpreter="latex");
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel,Interpreter="latex")  ; ylabel(CtrlVar.PlotsYaxisLabel,Interpreter="latex")
    colormap(othercolor("Mtemperaturemap",1028))
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
    fig=FindOrCreateFigure('C at the end of inversion') ; clf(fig)
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

    fig=FindOrCreateFigure('C at the beginning of inversion') ; clf(fig)
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

    fig=FindOrCreateFigure('Change in C during inversion run') ; clf(fig)
    PlotMeshScalarVariable(CtrlVar,MUA,log10(InvFinalValues.C)-log10(InvStartValues.C));
    title('log10(InvFinalValues.C)-log10(Cstart)') ;
    cbar=colorbar; title(cbar, '($\mathrm{m}\,\mathrm{yr}^{-1}\,\mathrm{kPa}^{-m}$)','interpreter','latex');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel,Interpreter="latex")  ; ylabel(CtrlVar.PlotsYaxisLabel,Interpreter="latex")
    colormap(othercolor("Mtemperaturemap",1028))
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
    title(cbar, '(m)');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    colormap(othercolor("Mdarkterrain",32))

    cbar=UaPlots(CtrlVar,MUA,F,InvStartValues.B,FigureTitle="B start");
    title('Bstart')
    title(cbar, '(m)')
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    colormap(othercolor("Mdarkterrain",32))

    cbar=UaPlots(CtrlVar,MUA,F,InvFinalValues.B-InvStartValues.B,FigureTitle="B final - B start");
    title('InvFinalValues.B-Bstart') ;
    title(cbar, '(m)');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    colormap(othercolor("Mdarkterrain",32))

    AspectRatio=1;
    fig=FindOrCreateFigure("sbB");  clf(fig)
    Plot_sbB(CtrlVar,MUA,[],[],F.B,[],[],AspectRatio) ; 
    title('B')
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    

end






%%

[~,~,tb] = CalcBasalTraction(CtrlVar,UserVar,MUA,F) ;
tb(tb<eps)=nan ;
cbar=UaPlots(CtrlVar,MUA,F,tb) ;
title('Basal drag, $\Vert \mathbf{t}_b \Vert$ ','interpreter','latex') ;
title(cbar, '($\mathrm{kPa}$)','interpreter','latex');
set(gca,'ColorScale','log')
clim([1 1000])
CM=cmocean('balanced') ; colormap(CM);
%%
% uAdjoint vAdjoint
if isprop(InvFinalValues,'uAdjoint')
    if ~isempty(InvFinalValues.uAdjoint)
        fig=FindOrCreateFigure('Adjoint variables') ; clf(fig)
        T=tiledlayout("flow");

        nexttile
        cbar=UaPlots(CtrlVar,MUA,F,InvFinalValues.uAdjoint,CreateNewFigure=false);
        title(' u Adjoint variable')
        CL=clim;
        if min(CL)< 0 && max(CL) > 0
            CM=cmocean('balanced',25,'pivot',0) ; colormap(CM);
        end

        nexttile
        cbar=UaPlots(CtrlVar,MUA,F,InvFinalValues.vAdjoint,CreateNewFigure=false);

        title(' v Adjoint variable')
        T.Padding="tight";   T.TileSpacing="tight";

        if min(CL)< 0 && max(CL) > 0
            CM=cmocean('balanced',25,'pivot',0) ; colormap(CM);
        end

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

fig=FindOrCreateFigure("Speed misfit") ; clf(fig)
speedMeas=sqrt(Meas.us.^2+Meas.vs.^2);
speedCalc=sqrt(F.ub.^2+F.vb.^2) ;
ErrSpeed=sqrt(usError.^2+vsError.^2);

T=tiledlayout("flow");

nexttile
cbar=UaPlots(CtrlVar,MUA,F,speedMeas,CreateNewFigure=false) ; title('Measured speed') ; set(gca,'ColorScale','log')
title(cbar,"$\|\mathbf{v}_\mathrm{Meas}\|$",interpreter="latex")
CL=clim; 

nexttile
cbar=UaPlots(CtrlVar,MUA,F,speedCalc,CreateNewFigure=false) ; title('Modelled speed') ; set(gca,'ColorScale','log')
title(cbar,"$\|\mathbf{v}_\mathrm{Modelled}\|$",interpreter="latex")
clim(CL);
subtitle("(Same colorbar scale as for measured speed)")

nexttile
cbar=UaPlots(CtrlVar,MUA,F,ErrSpeed,CreateNewFigure=false) ; title('Speed measurement error') ; set(gca,'ColorScale','log')
title(cbar,"error",interpreter="latex")


nexttile
D=speedMeas-speedCalc ;
cbar=UaPlots(CtrlVar,MUA,F,D,CreateNewFigure=false) ; title('Measured speed - modelled speed') ; set(gca,'ColorScale','lin')
title(cbar,"$\|\mathbf{v}_\mathrm{Meas}\|-\|\mathbf{v}_{\mathrm{Modelled}}\|$",interpreter="latex")
T.Padding="tight";   T.TileSpacing="tight";


%%
fig=FindOrCreateFigure('velocity misfit') ; clf(fig)


T=tiledlayout("flow");

nexttile
% Kplot=Kplot+1;     subplot(Iplot,Jplot,Kplot);
QuiverColorGHG(x,y,(us-Meas.us)./usError,(vs-Meas.vs)./vsError,CtrlVar);
title('((us-Meas.us)/usError,(vs-Meas.vs)/vsError)') ;
hold on
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
PlotMuaBoundary(CtrlVar,MUA,'b')  ;
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)

nexttile
%Kplot=Kplot+1;     subplot(Iplot,Jplot,Kplot);
QuiverColorGHG(x,y,us-Meas.us,vs-Meas.vs,CtrlVar); axis equal ; title('(us-Meas.us,v-Meas.vs)') ;
hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
PlotMuaBoundary(CtrlVar,MUA,'b')  ; xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)

if ~isempty(Meas.dhdt)  && contains(CtrlVar.Inverse.Measurements,"-dhdt")

    axdhdt=nexttile;
    [UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs);

    %Kplot=Kplot+1;
    %subplot(Iplot,Jplot,Kplot);
    PlotMeshScalarVariable(CtrlVar,MUA,(dhdt-Meas.dhdt)./dhdtError);
    hold on ;
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    PlotMuaBoundary(CtrlVar,MUA,'b')  ;
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)
    title('(dh/dt-Meas.dhdt)/dhdtError') ;
    CM=cmocean('balanced',25,'pivot',0) ; colormap(axdhdt,CM);

end

nexttile
%Kplot=Kplot+1;     subplot(Iplot,Jplot,Kplot);
[~,~,QuiverPar]=QuiverColorGHG(x,y,Meas.us,Meas.vs,CtrlVar); axis equal ;
title('(Meas.us,Meas.vs)') ;
hold on ;
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
PlotMuaBoundary(CtrlVar,MUA,'b')  ;
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)

nexttile
%Kplot=Kplot+1;     subplot(Iplot,Jplot,Kplot);
QuiverPar.QuiverSameVelocityScalingsAsBefore=1;
QuiverColorGHG(x,y,us,vs,QuiverPar); axis equal ; title('(us,vs)') ;
hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
PlotMuaBoundary(CtrlVar,MUA,'b')  ;
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)
QuiverPar.QuiverSameVelocityScalingsAsBefore=0;



if ~isempty(Meas.dhdt)  && contains(CtrlVar.Inverse.Measurements,"-dhdt")

    axdhdt=nexttile;
    %Kplot=Kplot+1; subplot(Iplot,Jplot,Kplot);
    PlotMeshScalarVariable(CtrlVar,MUA,dhdt);
    title('(dh/dt modelled)') ;
    hold on ;
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    PlotMuaBoundary(CtrlVar,MUA,'b')  ;
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)
    CM=cmocean('balanced',25,'pivot',0) ; colormap(axdhdt,CM);
end

T.Padding="tight";   T.TileSpacing="tight";

%%
fig=FindOrCreateFigure("Modelled velocities") ; clf(fig)
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

[~,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs);

UaPlots(CtrlVar,MUA,F,dhdt,FigureTitle="dh/dt modelled")
title('Modelled $dh/dt$ (assuming plug flow)','interpreter','latex') ;
CL=clim;
if CL(1) < 0 && CL(2)>0
    CM=cmocean('balanced',25,'pivot',0) ; colormap(CM);
else
    CM=cmocean('balanced',25) ;
end

%%  Prior

if isscalar(Priors.AGlen)
    Priors.AGlen=Priors.AGlen+zeros(MUA.Nnodes,1);
end

cbar=UaPlots(CtrlVar,MUA,F,Priors.AGlen,FigureTitle="APrior") ;
set(gca,'ColorScale','log')
title(cbar, '($\mathrm{yr}^{-1}\,\mathrm{kPa}^{-n}$)','interpreter','latex');
title("$A_{\mathrm{Prior}}$",Interpreter="latex")



cbar=UaPlots(CtrlVar,MUA,F,Priors.C,FigureTitle="CPrior") ;
set(gca,'ColorScale','log')
title(cbar, '($\mathrm{m}\,\mathrm{yr}^{-1}\,\mathrm{kPa}^{-m}$)','interpreter','latex');
title("$C_{\mathrm{Prior}}$",Interpreter="latex")





%%

if CtrlVar.Inverse.TestAdjoint.isTrue

   

    if ~isempty(InvFinalValues.dJdAGlenTest) | all(isnan(InvFinalValues.dJdAGlenTest))
        IA=find(~isnan(InvFinalValues.dJdAGlenTest)) ;
        fprintf('------------------------------------ AGlen gradients ---------------------------------------------------------------------\n')
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

        figAgrad=FindOrCreateFigure("A gradient test") ;  clf(figAgrad)
        plot(InvFinalValues.dJdAGlen,InvFinalValues.dJdAGlenTest,"or") ;
        hold on
        plot(InvFinalValues.dJdAGlen,InvFinalValues.dJdAGlen,"--k") ;
        axis equal tight ; 
        xlabel("Adjoint $dJ/dA$",Interpreter="latex")  ;
        ylabel("Finite difference $dJ/dA$",Interpreter="latex")
        ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
        axis on ; axis equal tight ; box off
        title("Comparision betweenadjoint and finite-differences gradient calculations")
        set(gcf,'Color','white')



    end

    if ~isempty(InvFinalValues.dJdCTest) | all(isnan(InvFinalValues.dJdCTest))

        IC=find(~isnan(InvFinalValues.dJdCTest)) ;

        fprintf('--------------------------------------- C gradients ----------------------------------------------------------------------\n')

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

        figCgrad=FindOrCreateFigure("C gradient test") ;  clf(figCgrad)
        plot(InvFinalValues.dJdC,InvFinalValues.dJdCTest,"or") ;
        hold on
        plot(InvFinalValues.dJdC,InvFinalValues.dJdC,"--k") ;
        axis equal tight  ; xlabel("Adjoint $dJ/dC$",Interpreter="latex")  ;
        ylabel("Finite difference $dJ/dC$",Interpreter="latex")
        ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
        axis on ; axis equal tight; 
        box off
        title("Comparision betweenadjoint and finite-differences gradient calculations")
        set(gcf,'Color','white')
    end

    if ~isempty(InvFinalValues.dJdBTest) | all(isnan(InvFinalValues.dJdBTest))
        fprintf('--------------------------------------- B gradients ----------------------------------------------------------------------\n')

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

        figBgrad=FindOrCreateFigure("B gradient test") ;  clf(figBgrad)
        plot(InvFinalValues.dJdB,InvFinalValues.dJdBTest,"or") ;
        hold on
        plot(InvFinalValues.dJdB,InvFinalValues.dJdB,"--k") ;
        axis equal ; xlabel("Adjoint $dJ/dB$",Interpreter="latex")  ;
        ylabel("Finite difference $dJ/dB$",Interpreter="latex")
        ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
        axis on ; axis equal tight ; box off
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

        IFigC=FindOrCreateFigure("dJ/dC gradient test over mesh") ; clf(IFigC);

        TileC=tiledlayout("flow");
        nexttile

        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdC) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('dJdC Adjoint gradient')

        nexttile
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdCTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('$dJ/dC$ Brute force gradient','interpreter','latex')


        nexttile
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdC-InvFinalValues.dJdCTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('Difference between adjoint and brute force derivatives')


        nexttile
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdC./InvFinalValues.dJdCTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('Ratio between adjoint and brute force derivatives')

        TileC.TileSpacing='tight';
        TileC.Padding='tight';
        %IFigC.Position=[948.43 41.571 1246.3 1115.4];
        %%
    end


    if ~(isempty(InvFinalValues.dJdAGlen) && isempty(InvFinalValues.dJdAGlenTest))


        IFigAGlen=FindOrCreateFigure("dJ/dA Test over mesh") ; clf(IFigAGlen);

        TileA=tiledlayout("flow");

        nexttile
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlen) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('dJdAGlen Adjoint gradient')

        nexttile
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlenTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('dJdAGlen Brute force gradient')


        nexttile
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlen-InvFinalValues.dJdAGlenTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('Difference between adjoint and brute force derivatives')

        nexttile
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlen./InvFinalValues.dJdAGlenTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('Ratio between adjoint and brute force derivatives')
        
        TileA.TileSpacing='tight';
        TileA.Padding='tight';
        % IFigAGlen.Position=[1.5714 41.571 1096 1115.4];
        %%
    end


    if ~(isempty(InvFinalValues.dJdB) && isempty(InvFinalValues.dJdBTest))

        %%
        IFigb=FindOrCreateFigure("dJ/dB gradient test over mesh") ; clf(IFigb)
        TileB=tiledlayout("flow");
        nexttile
        cbar=UaPlots(CtrlVar,MUA,F,InvFinalValues.dJdB,PlotUnderMesh=true,CreateNewFigure=false);
        title('$dJ/dB$ Adjoint gradient')

        nexttile
        cbar=UaPlots(CtrlVar,MUA,F,InvFinalValues.dJdBTest,PlotUnderMesh=true,CreateNewFigure=false);
        title('$dJ/dB$ Brute force gradient',Interpreter='latex')

        nexttile
        cbar=UaPlots(CtrlVar,MUA,F,InvFinalValues.dJdB-InvFinalValues.dJdBTest,PlotUnderMesh=true,CreateNewFigure=false);
        title('Difference between adjoint and brute force derivatives')

        nexttile
        UaPlots(CtrlVar,MUA,F,InvFinalValues.dJdB./InvFinalValues.dJdBTest,PlotUnderMesh=true,CreateNewFigure=false) ;
        title('Ratio between adjoint and brute force derivatives')

      %  IFigb.Position=[1.5714 41.571 1096 1115.4];
        TileB.TileSpacing='tight';
        TileB.Padding='tight';
        %%
    end

else

    if ~isempty(InvFinalValues.dJdAGlen)

        fig=FindOrCreateFigure('dJdA'); clf(fig) ;
        UaPlots(CtrlVar,MUA,F,InvFinalValues.dJdAGlen,CreateNewFigure=false);
        title('$dJ/dA$','interpreter','latex')
        cl=clim;
        if min(cl) <0 && max(cl)> 0
            CM=cmocean('balanced',25,'pivot',0) ; colormap(fig,CM);
        else
            CM=cmocean('balanced',25) ; colormap(fig,CM);
        end
    end

    if ~isempty(InvFinalValues.dJdC)

        fig=FindOrCreateFigure("dJdC"); clf(fig)
        UaPlots(CtrlVar,MUA,F,InvFinalValues.dJdC,CreateNewFigure=false);
        title('$dJ/dC$','interpreter','latex');
        cl=clim;
        if min(cl) <0 && max(cl)> 0
            CM=cmocean('balanced',25,'pivot',0) ; colormap(fig,CM);
        else
            CM=cmocean('balanced',25) ; colormap(fig,CM);
        end
    end

    if ~isempty(InvFinalValues.dJdB)

        fig=FindOrCreateFigure("dJdB"); clf(fig)
        UaPlots(CtrlVar,MUA,F,InvFinalValues.dJdB,CreateNewFigure=false);
        title('$dJ/dB$','interpreter','latex');
        cl=clim;
        if min(cl) <0 && max(cl)> 0
            CM=cmocean('balanced',25,'pivot',0) ; colormap(fig,CM);
        else
            CM=cmocean('balanced',25) ; colormap(fig,CM);
        end
    end


    %subplot(3,1,3) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dRdp) ; title('dRdp')


    %IFigGradients.Position=[1098.7 638.71 1096 518.29];
    %%
    CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;

    if contains(lower(CtrlVar.Inverse.InvertFor),'c')
        if ~isempty(Priors.TrueC)  && ~anynan(Priors.TrueC)


            figC=FindOrCreateFigure("True and estimated C"); clf(figC);

            T=tiledlayout("flow");

            nexttile
            UaPlots(CtrlVar,MUA,F,Priors.TrueC,CreateNewFigure=false) ;
            title("True C") ; set(gca,'ColorScale','log')

            nexttile
            UaPlots(CtrlVar,MUA,F,InvFinalValues.C,CreateNewFigure=false) ;
            title("Retrieved C") ; set(gca,'ColorScale','log')

            nexttile

            D=abs(Priors.TrueC-InvFinalValues.C) ;
            cbar=UaPlots(CtrlVar,MUA,F,D,CreateNewFigure=false) ;
            title('abs(True C - Retrieved C)') ; set(gca,'ColorScale','log')
            title(cbar,"$|C-\tilde{C}|$",interpreter="latex")

            nexttile
            UaPlots(CtrlVar,MUA,F,InvStartValues.C,CreateNewFigure=false);
            title("C at start of inversion") ; set(gca,'ColorScale','log')

            nexttile
            UaPlots(CtrlVar,MUA,F,Priors.C,CreateNewFigure=false) ;
            title('Prior C') ; set(gca,'ColorScale','log')


            nexttile
            UaPlots(CtrlVar,MUA,F,InvFinalValues.C-Priors.C,CreateNewFigure=false);
            title("Retrieved C -  Prior C ") ; set(gca,'ColorScale','log')

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

            nexttile
            UaPlots(CtrlVar,MUA,F,InvFinalValues.AGlen,CreateNewFigure=false) ;
            title('Retrieved AGlen') ; set(gca,'ColorScale','log')

            nexttile

            D=abs(Priors.TrueAGlen-InvFinalValues.AGlen) ;
            cbar=UaPlots(CtrlVar,MUA,F,D,CreateNewFigure=false) ;
            title('abs(True A -Retrieved A)') ; set(gca,'ColorScale','log')
            title(cbar,"$|A-\tilde{A}|$",interpreter="latex")

            nexttile
            UaPlots(CtrlVar,MUA,F,Priors.AGlen,CreateNewFigure=false) ;
            title('Prior AGlen') ; set(gca,'ColorScale','log')

            T.Padding="tight";   T.TileSpacing="tight";



        end
    end


    if contains(CtrlVar.Inverse.InvertFor,'-B-')

        if ~isempty(Priors.TrueB)  && ~anynan(Priors.TrueB)

            %% B

            PlotBedrockInversionFields(CtrlVar,MUA,F,Priors,InvFinalValues,InvStartValues,Meas)
           
            
          
        end

        if ~isempty(Priors.Trueh)  && ~anynan(Priors.Trueh)
            %% h

            figh=FindOrCreateFigure("True and estimated h"); clf(figh)
            TB=tiledlayout("flow") ;

            nexttile
            UaPlots(CtrlVar,MUA,F,Priors.Trueh,CreateNewFigure=false);
            title('True h')

            nexttile
            UaPlots(CtrlVar,MUA,F,F.h,CreateNewFigure=false);
            title('Retrieved h')

            nexttile
            UaPlots(CtrlVar,MUA,F,F.h-Priors.Trueh,CreateNewFigure=false);
            title('h estimated-true')

            nexttile
            [bStart,hStart]=Calc_bh_From_sBS(CtrlVar,MUA,F.s,InvStartValues.B,F.S,F.rho,F.rhow); %
            UaPlots(CtrlVar,MUA,F,hStart,CreateNewFigure=false);
            title("h at start of inversion")


            figB.Position=[500 200 900 800];
            TB.TileSpacing="tight";
            TB.Padding="tight";
            figh.Position=[500 200 900 800];

        end





        %%



    end







    %%
    if ~isempty(RunInfo.Inverse.J)

        fig=FindOrCreateFigure('Inverse Parameter Optimisation');
        clf(fig)
        hold off
        yyaxis left
        semilogy(RunInfo.Inverse.Iterations,RunInfo.Inverse.J,'-bo','LineWidth',2)
        ylabel('J','interpreter','latex')

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

            fig=FindOrCreateFigure('J=I+R');
            clf(fig)
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


end
