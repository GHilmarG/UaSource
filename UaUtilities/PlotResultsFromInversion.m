function PlotResultsFromInversion(UserVar,CtrlVar,MUA,BCs,F,~,~,InvStartValues,InvFinalValues,Priors,Meas,~,RunInfo)


%%
%
% PlotResultsFromInversion(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo)
%
% Does what it says on the tin.
%
%  Example:l
%
%  load InversionRestartFile
%  PlotResultsFromInversion(UserVarInRestartFile,CtrlVarInRestartFile,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);
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

fig=FindOrCreateFigure('Measuments') ; clf(fig)

Kplot=Kplot+1;    
subplot(Iplot,Jplot,Kplot)

PlotMeshScalarVariable(CtrlVar,MUA,Meas.us) ; 
hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');

xlabel(CtrlVar.PlotsXaxisLabel,'interpreter','latex');  
ylabel(CtrlVar.PlotsYaxisLabel,'interpreter','latex');
title('us Meas on numerical grid') ;

Kplot=Kplot+1;
subplot(Iplot,Jplot,Kplot)

PlotMeshScalarVariable(CtrlVar,MUA,Meas.vs) ; hold on ;
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');

xlabel(CtrlVar.PlotsXaxisLabel,'interpreter','latex');  
ylabel(CtrlVar.PlotsYaxisLabel,'interpreter','latex');
title('vs Meas on numerical grid') ;

if ~isempty(Meas.dhdt)  && contains(CtrlVar.Inverse.Measurements,"-dhdt")
    Kplot=Kplot+1;
    subplot(Iplot,Jplot,Kplot)
    PlotMeshScalarVariable(CtrlVar,MUA,Meas.dhdt) ; hold on ;
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    
    xlabel(CtrlVar.PlotsXaxisLabel,'interpreter','latex');
    ylabel(CtrlVar.PlotsYaxisLabel,'interpreter','latex');
    title('dh/dt Meas on numerical grid') ;
end

usError=sqrt(spdiags(Meas.usCov)); 
vsError=sqrt(spdiags(Meas.vsCov));
dhdtError=sqrt(spdiags(Meas.dhdtCov));


Kplot=Kplot+1;    
subplot(Iplot,Jplot,Kplot)


PlotMeshScalarVariable(CtrlVar,MUA,usError) ; hold on ;
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
title('us error on numerical grid') ;

Kplot=Kplot+1;
subplot(Iplot,Jplot,Kplot)
PlotMeshScalarVariable(CtrlVar,MUA,vsError) ; hold on ;
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
title('vs error on numerical grid') ;

if ~isempty(Meas.dhdt)  && contains(CtrlVar.Inverse.Measurements,"-dhdt")
    Kplot=Kplot+1;
    subplot(Iplot,Jplot,Kplot)
    PlotMeshScalarVariable(CtrlVar,MUA,dhdtError) ; hold on ;
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    title('dh/dt error on numerical grid') ;
    
end

    

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
    title('log10(InvFinalValues.AGlen)-log10(InvStartValues.AGlen)') ; cbar=colorbar; title(cbar, '($\mathrm{a}^{-1}$ $\mathrm{kPa}^{-3}$)',interpreter="latex");
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel,Interpreter="latex")  ; ylabel(CtrlVar.PlotsYaxisLabel,Interpreter="latex") 
    colormap(othercolor("Mtemperaturemap",1028))
    PlotMuaBoundary(CtrlVar,MUA,'k');
  

end

%%
if contains(upper(CtrlVar.Inverse.InvertFor),'C')
    
    ColorbarLimits=10.^[mean(log10(InvFinalValues.C))-4*std(log10(InvFinalValues.C))  mean(log10(InvFinalValues.C))+4*std(log10(InvFinalValues.C))];

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



if contains(CtrlVar.Inverse.InvertFor,'B')
    
    figure ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.B);
    title('InvFinalValues.B') ; cbar=colorbar; title(cbar, '(m)');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    
    figure ; PlotMeshScalarVariable(CtrlVar,MUA,InvStartValues.B);
    title('Bstart') ; cbar=colorbar; title(cbar, '(m)');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    
    figure ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.B-InvStartValues.B);
    title('InvFinalValues.B-Bstart') ; cbar=colorbar; title(cbar, '(m)');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    
    AspectRatio=1;
    figure ; Plot_sbB(CtrlVar,MUA,F.s,F.b,F.B,[],[],AspectRatio) ; title('F.s, F.b and F.B')
    
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
        subplot(1,2,1)
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.uAdjoint);
        hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title(' u Adjoint variable')
        
        subplot(1,2,2)
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.vAdjoint);
        hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title(' v Adjoint variable')
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

fig=FindOrCreateFigure('Speed misfit') ; clf(fig)
speedMeas=sqrt(Meas.us.^2+Meas.vs.^2);
speedCalc=sqrt(F.ub.^2+F.vb.^2) ;
ErrSpeed=sqrt(usError.^2+vsError.^2);


T=tiledlayout(2,2);

nexttile
cbar=UaPlots(CtrlVar,MUA,F,speedMeas,CreateNewFigure=false) ; title('Measured speed') ; set(gca,'ColorScale','log')
title(cbar,"$\|\mathbf{v}_\mathrm{Meas}\|$",interpreter="latex")

nexttile
cbar=UaPlots(CtrlVar,MUA,F,speedCalc,CreateNewFigure=false) ; title('Modelled speed') ; set(gca,'ColorScale','log')
title(cbar,"$\|\mathbf{v}_\mathrm{Modelled}\|$",interpreter="latex")

nexttile
cbar=UaPlots(CtrlVar,MUA,F,ErrSpeed,CreateNewFigure=false) ; title('Speed mesurement error') ; set(gca,'ColorScale','log')
title(cbar,"error",interpreter="latex")

nexttile
D=speedMeas-speedCalc ; 
cbar=UaPlots(CtrlVar,MUA,F,D,CreateNewFigure=false) ; title('Measured speed - modelled speed') ; set(gca,'ColorScale','log')
title(cbar,"$\|\mathbf{v}_\mathrm{Meas}\|-\|\mathbf{v}_{\mathrm{Modelled}}\|$",interpreter="latex")
T.Padding="tight";   T.TileSpacing="tight";


%%
fig=FindOrCreateFigure('velocity misfit') ; clf(fig)
Kplot=0;
T=tiledlayout;

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
    
    nexttile
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
     
    nexttile
    %Kplot=Kplot+1; subplot(Iplot,Jplot,Kplot);
    PlotMeshScalarVariable(CtrlVar,MUA,dhdt);
    title('(dh/dt modelled)') ;
    hold on ;
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    PlotMuaBoundary(CtrlVar,MUA,'b')  ;
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)
end

T.Padding="tight";   T.TileSpacing="tight";

%%
fig=FindOrCreateFigure('calculated velocities') ; clf(fig)
PlotBoundary(MUA.Boundary,MUA.connectivity,MUA.coordinates,CtrlVar,'k')
hold on
QuiverColorGHG(x,y,us,vs,QuiverPar); axis equal ; title("Calculated horizontal velocities") ;
hold on ;  
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,"r");
PlotCalvingFronts(CtrlVar,MUA,F,"b");
[~,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs); 

fig=FindOrCreateFigure('dh/dt calculated') ; clf(fig)
PlotBoundary(MUA.Boundary,MUA.connectivity,MUA.coordinates,CtrlVar,'k')
hold on
PlotMeshScalarVariable(CtrlVar,MUA,dhdt);
title('Calculated $dh/dt$ (assuming plug flow)','interpreter','latex') ;
hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
PlotCalvingFronts(CtrlVar,MUA,F,"b");

%%  Prior

if isscalar(Priors.AGlen)
Priors.AGlen=Priors.AGlen+zeros(MUA.Nnodes,1);
end

cbar=UaPlots(CtrlVar,MUA,F,Priors.AGlen,FigureTitle="log10(APrior)") ; 
set(gca,'ColorScale','log') 
title(cbar, '($\mathrm{yr}^{-1}\,\mathrm{kPa}^{-n}$)','interpreter','latex');
title("$A_{\mathrm{Prior}}$",Interpreter="latex")



cbar=UaPlots(CtrlVar,MUA,F,Priors.C,FigureTitle="log10(CPrior)") ; 
set(gca,'ColorScale','log') 
title(cbar, '($\mathrm{m}\,\mathrm{yr}^{-1}\,\mathrm{kPa}^{-m}$)','interpreter','latex');
title("$C_{\mathrm{Prior}}$",Interpreter="latex")





%%

if CtrlVar.Inverse.TestAdjoint.isTrue
    
    dJdp=InvFinalValues.dJdp;
    dJdpTest=InvFinalValues.dJdpTest;
%     fprintf('#Parameter  dJdp          dJdpTest      dJdp-dJdpTest     dJdp/dtdpTest \n')
%     
%     for ii=1:numel(iRange)
%         I=iRange(ii);
%         fprintf('%i %15g %15g  %15g  %15g \n',I,dJdp(I),dJdpTest(I),dJdp(I)-dJdpTest(I),dJdp(I)/dJdpTest(I))
%     end
    
  IA=find(~isnan(InvFinalValues.dJdAGlenTest)) ;
  fprintf('------------------------------------ AGlen gradients ---------------------------------------------------------------------\n')
    fprintf('#Node/Ele  dJdA          dJdATest      dJdA-dJdATest     dJdA/dtdATest \n')
    
    for ii=1:numel(IA)
        I=IA(ii);
        fprintf('%i %15g %15g  %15g  %15g \n',I,...
            InvFinalValues.dJdAGlen(I),...
            InvFinalValues.dJdAGlenTest(I),...
            InvFinalValues.dJdAGlen(I)-InvFinalValues.dJdAGlenTest(I),...
            InvFinalValues.dJdAGlen(I)/InvFinalValues.dJdAGlenTest(I))
    end
    
    
    IC=find(~isnan(InvFinalValues.dJdCTest)) ;
    
    fprintf('--------------------------------------- C gradients ----------------------------------------------------------------------\n')
    
    fprintf('#Node/Ele  dJdC          dJdCTest      dJdC-dJdCTest     dJdC/dtdCTest \n')
    
    for ii=1:numel(IC)
        I=IC(ii);
        fprintf('%i %15g %15g  %15g  %15g \n',I,...
            InvFinalValues.dJdC(I),...
            InvFinalValues.dJdCTest(I),...
            InvFinalValues.dJdC(I)-InvFinalValues.dJdCTest(I),...
            InvFinalValues.dJdC(I)/InvFinalValues.dJdCTest(I))
    end
    
    IB=find(~isnan(InvFinalValues.dJdBTest)) ;
    
    fprintf('--------------------------------------- B gradients ----------------------------------------------------------------------\n')
    
    fprintf('#Node/Ele  dJdB          dJdBTest      dJdB-dJdBTest     dJdB/dtdBTest \n')
    
    for ii=1:numel(IB)
        I=IB(ii);
        fprintf('%i %15g %15g  %15g  %15g \n',I,...
            InvFinalValues.dJdB(I),...
            InvFinalValues.dJdBTest(I),...
            InvFinalValues.dJdB(I)-InvFinalValues.dJdBTest(I),...
            InvFinalValues.dJdB(I)/InvFinalValues.dJdBTest(I))
    end
    
    
    fprintf('--------------------------------------------------------------------------------------------------------------------------\n')
    
    %[dJdp(iRange) dJdpTest(iRange)   dJdp(iRange)-dJdpTest(iRange) dJdp(iRange)./dJdpTest(iRange)]
    iRange=find(~isnan(dJdpTest));
    fprintf('Norm test: ||dJdpTest-dJdp||/||dJdp||= %g \n ',norm(dJdpTest(iRange)-dJdp(iRange))/norm(dJdp(iRange)))
    
    %%
    
    if ~(isempty(InvFinalValues.dJdC) && isempty(InvFinalValues.dJdCTest))
        
        IFigC=figure('Name','Inversion C','NumberTitle','off');
        
        
        
        subplot(2,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdC) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('dJdC Adjoint gradient')
        
        subplot(2,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdCTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('$dJ/dC$ Brute force gradient','interpreter','latex')
        
        
        subplot(2,2,3) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdC-InvFinalValues.dJdCTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('Difference between adjoint and brute force derivatives')
        
        subplot(2,2,4) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdC./InvFinalValues.dJdCTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('Ratio between adjoint and brute force derivatives')
        
        IFigC.Position=[948.43 41.571 1246.3 1115.4];
        %%
    end
    
    
    if ~(isempty(InvFinalValues.dJdAGlen) && isempty(InvFinalValues.dJdAGlenTest))
        
        IFigAGlen=figure('Name','Inversion AGlen','NumberTitle','off');
        
        
        
        subplot(2,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlen) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('dJdAGlen Adjoint gradient')
        
        subplot(2,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlenTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('dJdAGlen Brute force gradient')
        
        
        subplot(2,2,3) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlen-InvFinalValues.dJdAGlenTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('Difference between adjoint and brute force derivatives')
        
        subplot(2,2,4) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlen./InvFinalValues.dJdAGlenTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('Ratio between adjoint and brute force derivatives')
        
        IFigAGlen.Position=[1.5714 41.571 1096 1115.4];
        %%
    end
    
       
    if ~(isempty(InvFinalValues.dJdB) && isempty(InvFinalValues.dJdBTest))
        
        IFigAGlen=figure('Name','Inversion b','NumberTitle','off');
        
        
        
        subplot(2,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdB) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('dJdB Adjoint gradient')
        
        subplot(2,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdBTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('dJdB Brute force gradient')
        
        
        subplot(2,2,3) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdB-InvFinalValues.dJdBTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('Difference between adjoint and brute force derivatives')
        
        subplot(2,2,4) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdB./InvFinalValues.dJdBTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('Ratio between adjoint and brute force derivatives')
        
        IFigAGlen.Position=[1.5714 41.571 1096 1115.4];
        %%
    end
    
    
    
       
    if ~(isempty(InvFinalValues.dJdB) && isempty(InvFinalValues.dJdBTest))
        
        IFigAGlen=figure('Name','Inversion B','NumberTitle','off');
        
        
        
        subplot(2,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdB) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('dJdB Adjoint gradient')
        
        subplot(2,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdBTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('dJdB Brute force gradient')
        
        
        subplot(2,2,3) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdB-InvFinalValues.dJdBTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('Difference between adjoint and brute force derivatives')
        
        subplot(2,2,4) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdB./InvFinalValues.dJdBTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
        title('Ratio between adjoint and brute force derivatives')
        
        IFigAGlen.Position=[1.5714 41.571 1096 1115.4];
        %%
    end
    
    
    %%
else
    
    if ~isempty(InvFinalValues.dJdAGlen)
        
        fig=FindOrCreateFigure('dJdA'); clf(fig) ;
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlen) ; 
        title('$dJ/dA$','interpreter','latex')
        hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    end
    
    if ~isempty(InvFinalValues.dJdC)
        
        fig=FindOrCreateFigure('dJdC'); clf(fig)
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdC) ; 
        title('$dJ/dC$','interpreter','latex');
        hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    end
    
    if ~isempty(InvFinalValues.dJdB)
        
        fig=FindOrCreateFigure('dJdB'); clf(fig)
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdB) ;
        title('$dJ/dC$','interpreter','latex');
        hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
    end
    
    
    %subplot(3,1,3) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dRdp) ; title('dRdp')
    
    
    %IFigGradients.Position=[1098.7 638.71 1096 518.29];
    %%
    CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
    
    if contains(lower(CtrlVar.Inverse.InvertFor),'c')
        if ~isempty(Priors.TrueC)


            tFig1=FindOrCreateFigure("True and estimated C"); clf(tFig1 ); 
            
            T=tiledlayout(2,2);
            
            nexttile
            UaPlots(CtrlVar,MUA,F,Priors.TrueC,CreateNewFigure=false) ; 
            title('True C') ; set(gca,'ColorScale','log')

            nexttile
            UaPlots(CtrlVar,MUA,F,InvFinalValues.C,CreateNewFigure=false) ; 
            title('Retrieved C') ; set(gca,'ColorScale','log')

            nexttile
            
            D=abs(Priors.TrueC-InvFinalValues.C) ; 
            cbar=UaPlots(CtrlVar,MUA,F,D,CreateNewFigure=false) ; 
            title('abs(True C - Retrieved C)') ; set(gca,'ColorScale','log')
            title(cbar,"$|C-\tilde{C}|$",interpreter="latex")

            nexttile
            UaPlots(CtrlVar,MUA,F,Priors.C,CreateNewFigure=false) ; 
            title('Prior C') ; set(gca,'ColorScale','log')

            T.Padding="tight";   T.TileSpacing="tight";


        end
    end
    
    if contains(lower(CtrlVar.Inverse.InvertFor),'aglen')
        
        if ~isempty(Priors.TrueAGlen)



            tFig1=FindOrCreateFigure("True and estimated AGlen"); clf(tFig1 ); 
            
            T=tiledlayout(2,2);
            
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
        
        if ~isempty(Priors.TrueB)
            
            tFig1=figure('Name','True and estimated B','NumberTitle','off');
            subplot(1,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,Priors.TrueB) ; title('True B')
            hold on ; PlotMuaMesh(CtrlVar,MUA,[],'w') ;
            subplot(1,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.B) ; title('Retrieved B')
            hold on ; PlotMuaMesh(CtrlVar,MUA,[],'w') ; 
            tFig1.Units='normalized';
            tFig1.Position=[0.5 0.51 0.5 0.4];
            
            
            tFig2=figure('Name','Difference between true and estimated b','NumberTitle','off');
            
            
                  
            subplot(2,2,1)
            PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.B);
            %hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
            SetLabels(CtrlVar,"km","km","m");
            title("Inv. start field: "+CtrlVar.Inverse.InvertFor)
            colorbar('southoutside')
            
            subplot(2,2,2)
            PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.B);
            %hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
            SetLabels(CtrlVar,"km","km","m");
            title("Inverted: "+CtrlVar.Inverse.InvertFor)
            colorbar('southoutside')
            
            subplot(2,2,3)
            PlotMeshScalarVariable(CtrlVar,MUA,Priors.TrueB);            
            %hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
            SetLabels(CtrlVar,"km","km","m");
            title("True: "+CtrlVar.Inverse.InvertFor)
            colorbar('southoutside')
            
            subplot(2,2,4)
            PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.B-Priors.TrueB);
            %hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
            SetLabels(CtrlVar,"km","km","m");
            title("Estimated-True: "+CtrlVar.Inverse.InvertFor)
            colorbar('southoutside')
            
            tFig2.Units='normalized';
            tFig2.Position=[0.3 0.2 0.8 0.5];
            
            tFigTh=figure('Name','Difference between true and estimated h','NumberTitle','off');
            
            subplot(1,3,1)
            PlotMeshScalarVariable(CtrlVar,MUA,F.s-InvFinalValues.B);
            hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
            SetLabels(CtrlVar,"km","km","m");
            title("Inverted: h")
            colorbar('off')
            colorbar('southoutside')
            
            subplot(1,3,2)
            PlotMeshScalarVariable(CtrlVar,MUA,F.s-Priors.TrueB)  ; % this may need to be adjusted
            hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
            SetLabels(CtrlVar,"km","km","m");
            title("True: h")
            colorbar('off')
            colorbar('southoutside')
            
            subplot(1,3,3)
            PlotMeshScalarVariable(CtrlVar,MUA,(F.s-InvFinalValues.B)-(F.s-Priors.TrueB));
            hold on ; PlotGroundingLines(CtrlVar,MUA,F.GF,GLgeo,xGL,yGL,'r');
           SetLabels(CtrlVar,"km","km","m");
            title("Inverted-True: h ")
            colorbar('off')
            colorbar('southoutside')
            
            tFigTh.Units='normalized';
            
            tFigTh.Position=[0.1 0.4 0.8 0.5];
            
        end
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
%%
fprintf('done.\n')

end
