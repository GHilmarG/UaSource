function PlotResultsFromInversion(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo)


%%
%
% PlotResultsFromInversion(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo)
%
% Does what it says on the tin.
%
%  Example:
%
%  load InversionRestartFile
%  PlotResultsFromInversion(UserVarInRestartFile,CtrlVarInRestartFile,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);
%
%%

fprintf(' Plotting results from inversion...')

CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
us=F.ub+F.ud; vs=F.vb+F.vd;

x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);

GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar); xGL=[] ; yGL=[] ;

%%

if ~isempty(Meas.dhdt)
     Iplot=2 ; Jplot=3;
else
     Iplot=2 ; Jplot=2;
end
Kplot=0;    

fig=FindOrCreateFigure('Measuments') ;

Kplot=Kplot+1;    
subplot(Iplot,Jplot,Kplot)

PlotMeshScalarVariable(CtrlVar,MUA,Meas.us) ; 
hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');

xlabel(CtrlVar.PlotsXaxisLabel,'interpreter','latex');  
ylabel(CtrlVar.PlotsYaxisLabel,'interpreter','latex');
title('us Meas on numerical grid') ;

Kplot=Kplot+1;
subplot(Iplot,Jplot,Kplot)

PlotMeshScalarVariable(CtrlVar,MUA,Meas.vs) ; hold on ;
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');

xlabel(CtrlVar.PlotsXaxisLabel,'interpreter','latex');  
ylabel(CtrlVar.PlotsYaxisLabel,'interpreter','latex');
title('vs Meas on numerical grid') ;

if ~isempty(Meas.dhdt)
    Kplot=Kplot+1;
    subplot(Iplot,Jplot,Kplot)
    PlotMeshScalarVariable(CtrlVar,MUA,Meas.dhdt) ; hold on ;
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    
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
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
title('us error on numerical grid') ;

Kplot=Kplot+1;
subplot(Iplot,Jplot,Kplot)
PlotMeshScalarVariable(CtrlVar,MUA,vsError) ; hold on ;
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
title('vs error on numerical grid') ;

if ~isempty(Meas.dhdt)
    Kplot=Kplot+1;
    subplot(Iplot,Jplot,Kplot)
    PlotMeshScalarVariable(CtrlVar,MUA,dhdtError) ; hold on ;
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    title('dh/dt error on numerical grid') ;
    
end

    

%%
if contains(upper(CtrlVar.Inverse.InvertFor),'A')
    
    fig=FindOrCreateFigure('A at the end of inversion') ;
    PlotMeshScalarVariable(CtrlVar,MUA,log10(InvFinalValues.AGlen));
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    CtrlVar.PlotNodes=0 ; % PlotMuaMesh(CtrlVar,MUA,[],'k') ; 
    title('log10(InvFinalValues.AGlen)') ; cbar=colorbar; title(cbar, '($\mathrm{a}^{-1}$ $\mathrm{kPa}^{-3}$)',interpreter="latex");
    
    fig=FindOrCreateFigure('A at the beginning of inversion') ;
    PlotMeshScalarVariable(CtrlVar,MUA,log10(InvStartValues.AGlen));
    title('log10(Astart)') ; cbar=colorbar; title(cbar, '($\mathrm{a}^{-1}$ $\mathrm{kPa}^{-3}$)',interpreter="latex");
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    
    fig=FindOrCreateFigure('Change in A during inversion run') ;
    PlotMeshScalarVariable(CtrlVar,MUA,log10(InvFinalValues.AGlen)-log10(InvStartValues.AGlen));
    title('log10(InvFinalValues.AGlen)-log10(InvStartValues.AGlen)') ; cbar=colorbar; title(cbar, '($\mathrm{a}^{-1}$ $\mathrm{kPa}^{-3}$)',interpreter="latex");
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
end

%%
if contains(upper(CtrlVar.Inverse.InvertFor),'C')
    
   
    fig=FindOrCreateFigure('C at the end of inversion') ;
    PlotMeshScalarVariable(CtrlVar,MUA,log10(InvFinalValues.C));
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    CtrlVar.PlotNodes=0 ; % PlotMuaMesh(CtrlVar,MUA,[],'k') ; 
    title('log10(InvFinalValues.C)','interpreter','latex');
    cbar=colorbar; title(cbar, '($m\,\mathrm{yr}^{-1}\,\mathrm{kPa}^{-m}$)','interpreter','latex');
    
    fig=FindOrCreateFigure('C at the beginning of inversion') ;
    PlotMeshScalarVariable(CtrlVar,MUA,log10(InvStartValues.C));
    title('log10(Cstart)') ; 
    cbar=colorbar; title(cbar, '($m\,\mathrm{yr}^{-1}\,\mathrm{kPa}^{-m}$)','interpreter','latex');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    
    fig=FindOrCreateFigure('Change in C during inversion run') ;
    PlotMeshScalarVariable(CtrlVar,MUA,log10(InvFinalValues.C)-log10(InvStartValues.C));
    title('log10(InvFinalValues.C)-log10(Cstart)') ; 
    cbar=colorbar; title(cbar, '($m\,\mathrm{yr}^{-1}\,\mathrm{kPa}^{-m}$)','interpreter','latex');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
end

if contains(CtrlVar.Inverse.InvertFor,'b')
    
    figure ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.b);
    title('InvFinalValues.b') ; cbar=colorbar; title(cbar, '(m)');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    
    figure ; PlotMeshScalarVariable(CtrlVar,MUA,InvStartValues.b);
    title('bstart') ; cbar=colorbar; title(cbar, '(m)');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    
    figure ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.b-InvStartValues.b);
    title('InvFinalValues.b-bstart') ; cbar=colorbar; title(cbar, '(m)');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    
    %[TRI,DT,LightHandle]=Plot_sbB(CtrlVar,MUA,s,b,B,TRI,DT,AspectRatio,ViewAndLight,LightHandle,sCol,bCol,BCol);
    AspectRatio=1;
    figure ; Plot_sbB(CtrlVar,MUA,F.s,F.b,F.B,[],[],AspectRatio) ; title('F.s, F.b and F.B')
end



if contains(CtrlVar.Inverse.InvertFor,'B')
    
    figure ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.B);
    title('InvFinalValues.B') ; cbar=colorbar; title(cbar, '(m)');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    
    figure ; PlotMeshScalarVariable(CtrlVar,MUA,InvStartValues.B);
    title('Bstart') ; cbar=colorbar; title(cbar, '(m)');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    
    figure ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.B-InvStartValues.B);
    title('InvFinalValues.B-Bstart') ; cbar=colorbar; title(cbar, '(m)');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    
    AspectRatio=1;
    figure ; Plot_sbB(CtrlVar,MUA,F.s,F.b,F.B,[],[],AspectRatio) ; title('F.s, F.b and F.B')
    
end






%%
F.C=InvFinalValues.C; % this should already have been updated internally in Ua
[tbx,tby,tb] = CalcBasalTraction(CtrlVar,UserVar,MUA,F) ;
fig=FindOrCreateFigure('Basal traction') ;
PlotMeshScalarVariable(CtrlVar,MUA,tb) ;
title('Basal drag, $\Vert \mathbf{t}_b \Vert$ ','interpreter','latex') ;
cbar=colorbar; title(cbar, '($\mathrm{kPa}$)','interpreter','latex');
hold on
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');

%%
% uAdjoint vAdjoint
if isprop(InvFinalValues,'uAdjoint')
    if ~isempty(InvFinalValues.uAdjoint)
        fig=FindOrCreateFigure('Adjoint variables') ;
        subplot(1,2,1)
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.uAdjoint);
        hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
        title(' u Adjoint variable')
        
        subplot(1,2,2)
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.vAdjoint);
        hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
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
    wsError=sqrt(spdiags(Meas.wsCov));
end
if ~exist('GLgeo','var')
    GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar); xGL=[] ; yGL=[] ;
end

%%

fig=FindOrCreateFigure('Speed misfit') ;
speedMeas=sqrt(Meas.us.^2+Meas.vs.^2);
speedCalc=sqrt(F.ub.^2+F.vb.^2) ;
ErrSpeed=sqrt(usError.^2+vsError.^2); 

subplot(2,2,1)
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,log10(speedMeas)) ; title('log10(measured speed)') 
hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
PlotMuaBoundary(CtrlVar,MUA,'b')  ; xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
cAxisMeas=caxis; 

subplot(2,2,2)
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,log10(speedCalc)) ; title('log10(calculated speed)') 
hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
PlotMuaBoundary(CtrlVar,MUA,'b')  ; xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
caxis(cAxisMeas);

subplot(2,2,3)
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,log10(ErrSpeed)) ; title('log10(Meas error in speed)') 
hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
PlotMuaBoundary(CtrlVar,MUA,'b')  ; xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);

subplot(2,2,4)
[~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,(speedMeas-speedCalc)./ErrSpeed) ; title('speed residuals: (meas-calc)/error') 
hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
PlotMuaBoundary(CtrlVar,MUA,'b')  ; xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
%%
fig=FindOrCreateFigure('velocity misfit') ;
Kplot=0;

Kplot=Kplot+1;    
subplot(Iplot,Jplot,Kplot);
QuiverColorGHG(x,y,(us-Meas.us)./usError,(vs-Meas.vs)./vsError,CtrlVar);
title('((us-Meas.us)/usError,(vs-Meas.vs)/vsError)') ;
hold on
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
PlotMuaBoundary(CtrlVar,MUA,'b')  ;
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)

Kplot=Kplot+1;    
subplot(Iplot,Jplot,Kplot);
QuiverColorGHG(x,y,us-Meas.us,vs-Meas.vs,CtrlVar); axis equal ; title('(us-Meas.us,v-Meas.vs)') ;
hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
PlotMuaBoundary(CtrlVar,MUA,'b')  ; xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)

if ~isempty(Meas.dhdt)
    
    [UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs);
     
    Kplot=Kplot+1;
    subplot(Iplot,Jplot,Kplot);
    PlotMeshScalarVariable(CtrlVar,MUA,(dhdt-Meas.dhdt)./dhdtError);
    hold on ;
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    PlotMuaBoundary(CtrlVar,MUA,'b')  ;
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)
    title('(dh/dt-Meas.dhdt)/dhdtError') ;
    
end


Kplot=Kplot+1;    
subplot(Iplot,Jplot,Kplot);
[~,~,QuiverPar]=QuiverColorGHG(x,y,Meas.us,Meas.vs,CtrlVar); axis equal ;
title('(Meas.us,Meas.vs)') ;
hold on ;
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
PlotMuaBoundary(CtrlVar,MUA,'b')  ;
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)

Kplot=Kplot+1;    
subplot(Iplot,Jplot,Kplot);
QuiverPar.QuiverSameVelocityScalingsAsBefore=1;
QuiverColorGHG(x,y,us,vs,QuiverPar); axis equal ; title('(us,vs)') ;
hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
PlotMuaBoundary(CtrlVar,MUA,'b')  ;
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)
QuiverPar.QuiverSameVelocityScalingsAsBefore=0;


if ~isempty(Meas.dhdt)
     
    Kplot=Kplot+1;
    subplot(Iplot,Jplot,Kplot);
    PlotMeshScalarVariable(CtrlVar,MUA,dhdt);
    title('(dh/dt modelled)') ;
    hold on ;
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    PlotMuaBoundary(CtrlVar,MUA,'b')  ;
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)
end



%%
fig=FindOrCreateFigure('calculated velocities') ;
PlotBoundary(MUA.Boundary,MUA.connectivity,MUA.coordinates,CtrlVar,'k')
hold on
QuiverColorGHG(x,y,us,vs,QuiverPar); axis equal ; title('Calculated horizontal velocities') ;
hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');

[UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs); 

fig=FindOrCreateFigure('dh/dt calculated') ;
PlotBoundary(MUA.Boundary,MUA.connectivity,MUA.coordinates,CtrlVar,'k')
hold on
PlotMeshScalarVariable(CtrlVar,MUA,dhdt);
title('Calculated $dh/dt$ (assuming plug flow)','interpreter','latex') ;
hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');

%%  % Difference in speed

SpeedMeas=sqrt(Meas.us.^2+Meas.vs.^2);
SpeedCalc=sqrt(us.^2+vs.^2);

SpeedDiff=100*(SpeedCalc-SpeedMeas)./SpeedMeas;
fig=FindOrCreateFigure('Normalized speed misfit') ;

PlotMeshScalarVariable(CtrlVar,MUA,SpeedDiff);
hold on 
hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
title('100*(SpeedCalc-SpeedMeas)./SpeedMeas')
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)

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
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
        title('dJdC Adjoint gradient')
        
        subplot(2,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdCTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
        title('$dJ/dC$ Brute force gradient','interpreter','latex')
        
        
        subplot(2,2,3) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdC-InvFinalValues.dJdCTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
        title('Difference between adjoint and brute force derivatives')
        
        subplot(2,2,4) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdC./InvFinalValues.dJdCTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
        title('Ratio between adjoint and brute force derivatives')
        
        IFigC.Position=[948.43 41.571 1246.3 1115.4];
        %%
    end
    
    
    if ~(isempty(InvFinalValues.dJdAGlen) && isempty(InvFinalValues.dJdAGlenTest))
        
        IFigAGlen=figure('Name','Inversion AGlen','NumberTitle','off');
        
        
        
        subplot(2,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlen) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
        title('dJdAGlen Adjoint gradient')
        
        subplot(2,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlenTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
        title('dJdAGlen Brute force gradient')
        
        
        subplot(2,2,3) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlen-InvFinalValues.dJdAGlenTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
        title('Difference between adjoint and brute force derivatives')
        
        subplot(2,2,4) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlen./InvFinalValues.dJdAGlenTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
        title('Ratio between adjoint and brute force derivatives')
        
        IFigAGlen.Position=[1.5714 41.571 1096 1115.4];
        %%
    end
    
       
    if ~(isempty(InvFinalValues.dJdB) && isempty(InvFinalValues.dJdBTest))
        
        IFigAGlen=figure('Name','Inversion b','NumberTitle','off');
        
        
        
        subplot(2,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdB) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
        title('dJdB Adjoint gradient')
        
        subplot(2,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdBTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
        title('dJdB Brute force gradient')
        
        
        subplot(2,2,3) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdB-InvFinalValues.dJdBTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
        title('Difference between adjoint and brute force derivatives')
        
        subplot(2,2,4) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdB./InvFinalValues.dJdBTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
        title('Ratio between adjoint and brute force derivatives')
        
        IFigAGlen.Position=[1.5714 41.571 1096 1115.4];
        %%
    end
    
    
    
       
    if ~(isempty(InvFinalValues.dJdB) && isempty(InvFinalValues.dJdBTest))
        
        IFigAGlen=figure('Name','Inversion B','NumberTitle','off');
        
        
        
        subplot(2,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdB) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
        title('dJdB Adjoint gradient')
        
        subplot(2,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdBTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
        title('dJdB Brute force gradient')
        
        
        subplot(2,2,3) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdB-InvFinalValues.dJdBTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
        title('Difference between adjoint and brute force derivatives')
        
        subplot(2,2,4) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdB./InvFinalValues.dJdBTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        hold on ;  [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
        title('Ratio between adjoint and brute force derivatives')
        
        IFigAGlen.Position=[1.5714 41.571 1096 1115.4];
        %%
    end
    
    
    %%
else
    
    if ~isempty(InvFinalValues.dJdAGlen)
        
        fig=FindOrCreateFigure('dJdA');
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlen) ; 
        title('$dJ/dA$','interpreter','latex')
        hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    end
    
    if ~isempty(InvFinalValues.dJdC)
        
        fig=FindOrCreateFigure('dJdC');
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdC) ; 
        title('$dJ/dC$','interpreter','latex');
        hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    end
    
    if ~isempty(InvFinalValues.dJdB)
        
        fig=FindOrCreateFigure('dJdB');
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdB) ;
        title('$dJ/dC$','interpreter','latex');
        hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    end
    
    
    %subplot(3,1,3) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dRdp) ; title('dRdp')
    
    
    %IFigGradients.Position=[1098.7 638.71 1096 518.29];
    %%
    CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
    
    if contains(lower(CtrlVar.Inverse.InvertFor),'c')
        if ~isempty(Priors.TrueC)
            tFig1=figure('Name','True and estimated C','NumberTitle','off');

            subplot(1,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,Priors.TrueC) ; title('True C')
            hold on ; PlotMuaMesh(CtrlVar,MUA,[],'w');
            subplot(1,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.C) ; title('Retrieved C')
            hold on ; PlotMuaMesh(CtrlVar,MUA,[],'w');
            tFig1.Units='normalized';
            tFig1.Position=[0.5 0.5 0.5 0.4];
            
            if  ~CtrlVar.CisElementBased
                tFig2=figure('Name','Difference between true and estimated','NumberTitle','off');
                %PlotMeshScalarVariable(CtrlVar,MUA,Priors.TrueC-InvFinalValues.C);
                
                subplot(1,3,1)
                [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.C); 
                xlabel('x') ; ylabel('y') ; title('Inverted slipperiness')
                colorbar('southoutside')
                
                subplot(1,3,2)
                [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,Priors.TrueC); 
                xlabel('x') ; ylabel('y') ; title('True slipperiness')
                colorbar('southoutside')
                
                subplot(1,3,3)
                [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.C-Priors.TrueC);
                xlabel('x') ; ylabel('y') ; title('True slipperiness')
                title('Slipperiness: True-Estimated')
                colorbar('southoutside')
                
                tFig2.Units='normalized';
                tFig2.Position=[0.1 0.2 0.8 0.5];
            end
        end
    end
    
    if contains(lower(CtrlVar.Inverse.InvertFor),'aglen')
        
        if ~isempty(Priors.TrueAGlen)
            tFig1=figure('Name','True and estimated AGlen','NumberTitle','off');
            subplot(1,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,Priors.TrueAGlen) ; title('True AGlen')
            hold on ; PlotMuaMesh(CtrlVar,MUA,[],'w');
            subplot(1,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.AGlen) ; title('Retrieved AGlen')
            hold on ; PlotMuaMesh(CtrlVar,MUA,[],'w');
            tFig1.Units='normalized';
            tFig1.Position=[0.5 0.5 0.5 0.4];
            
            if  ~CtrlVar.AGlenisElementBased
                tFig2=figure('Name','Difference between true and estimated','NumberTitle','off');
                %PlotMeshScalarVariable(CtrlVar,MUA,Priors.TrueC-InvFinalValues.C);
                
                subplot(1,3,1)
                [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.AGlen);
                
                SetLabels("km","km","m");
                title('Inverted AGlen')
                colorbar('southoutside')
                
                subplot(1,3,2)
                [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,Priors.TrueAGlen);
                SetLabels("km","km","m");
                title('True AGlen')
                colorbar('southoutside')
                
                subplot(1,3,3)
                [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.AGlen-Priors.TrueAGlen);
                SetLabels("km","km","m");
                title('True AGlen')
                title('AGlen: True-Estimated')
                colorbar('southoutside')
                
                tFig2.Units='normalized';
                tFig2.Position=[0.2 0.2 0.8 0.5];
            end
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
            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.B);
            %hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
            SetLabels(CtrlVar,"km","km","m");
            title("Inv. start field: "+CtrlVar.Inverse.InvertFor)
            colorbar('southoutside')
            
            subplot(2,2,2)
            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.B);
            %hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
            SetLabels(CtrlVar,"km","km","m");
            title("Inverted: "+CtrlVar.Inverse.InvertFor)
            colorbar('southoutside')
            
            subplot(2,2,3)
            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,Priors.TrueB);            
            %hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
            SetLabels(CtrlVar,"km","km","m");
            title("True: "+CtrlVar.Inverse.InvertFor)
            colorbar('southoutside')
            
            subplot(2,2,4)
            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.B-Priors.TrueB);
            %hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
            SetLabels(CtrlVar,"km","km","m");
            title("Estimated-True: "+CtrlVar.Inverse.InvertFor)
            colorbar('southoutside')
            
            tFig2.Units='normalized';
            tFig2.Position=[0.3 0.2 0.8 0.5];
            
            tFigTh=figure('Name','Difference between true and estimated h','NumberTitle','off');
            
            subplot(1,3,1)
            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,F.s-InvFinalValues.B);
            hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
            SetLabels(CtrlVar,"km","km","m");
            title("Inverted: h")
            colorbar('off')
            colorbar('southoutside')
            
            subplot(1,3,2)
            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,F.s-Priors.TrueB)  ; % this may need to be adjusted
            hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
            SetLabels(CtrlVar,"km","km","m");
            title("True: h")
            colorbar('off')
            colorbar('southoutside')
            
            subplot(1,3,3)
            [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,(F.s-InvFinalValues.B)-(F.s-Priors.TrueB));
            hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
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
            ylabel('J','interpreter','latex')
            
            hold on
            semilogy(RunInfo.Inverse.Iterations,RunInfo.Inverse.I,'-gx')
            ylabel('J and I')
         
            yyaxis right
            semilogy(RunInfo.Inverse.Iterations,RunInfo.Inverse.R,'-r+')
            ylabel('R','interpreter','latex')
            legend('Objective function','I','R','Location','southwest','interpreter','latex')

        end
        
        yyaxis left
        xlabel('Inverse iteration','interpreter','latex') ;
        hold off
    end
    
    
end
%%
fprintf('done.\n')

end
