function AdjointResultsPlots(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo)


%%
CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
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
if contains(upper(CtrlVar.Inverse.InvertFor),'A')
    
    figure ; PlotMeshScalarVariable(CtrlVar,MUA,log10(InvFinalValues.AGlen));
    title('log10(InvFinalValues.AGlen)') ; cbar=colorbar; title(cbar, '(a^{-1} kPa^{-3})');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    
    figure ; PlotMeshScalarVariable(CtrlVar,MUA,log10(InvStartValues.AGlen));
    title('log10(Astart)') ; cbar=colorbar; title(cbar, '(a^{-1} kPa^{-3})');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    
    figure ; PlotMeshScalarVariable(CtrlVar,MUA,log10(InvFinalValues.AGlen)-log10(InvStartValues.AGlen));
    title('log10(InvFinalValues.AGlen)-log10(InvStartValues.AGlen)') ; cbar=colorbar; title(cbar, '(a^{-1} kPa^{-3})');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
end

%%
if contains(upper(CtrlVar.Inverse.InvertFor),'C')
    
    figure ; PlotMeshScalarVariable(CtrlVar,MUA,log10(InvFinalValues.C));
    title('log10(InvFinalValues.C)') ; cbar=colorbar; title(cbar, '(m/a/kPa^m)');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    
    figure ; PlotMeshScalarVariable(CtrlVar,MUA,log10(InvStartValues.C));
    title('log10(Cstart)') ; cbar=colorbar; title(cbar, '(m/a/kPa^m)');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
    
    figure ; PlotMeshScalarVariable(CtrlVar,MUA,log10(InvFinalValues.C)-log10(InvStartValues.C));
    title('log10(InvFinalValues.C)-log10(Cstart)') ; cbar=colorbar; title(cbar, '(m/a/kPa^m)');
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
end


%%
[tbx,tby,tb,beta2] = CalcBasalTraction(CtrlVar,MUA,F.ub,F.vb,InvFinalValues.C,InvFinalValues.m,GF) ;
figure
PlotMeshScalarVariable(CtrlVar,MUA,tb) ;
title(' tb ') ; cbar=colorbar; title(cbar, '(kPa)');
hold on
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');

%%
% uAdjoint vAdjoint
if isprop(InvFinalValues,'uAdjoint')
    if ~isempty(InvFinalValues.uAdjoint)
        figure ;
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
figure
subplot(2,2,1);
QuiverColorGHG(x,y,(us-Meas.us)./usError,(vs-Meas.vs)./vsError,CtrlVar);
title('((us-Meas.us)/usError,(vs-Meas.vs)/vsError)') ;
hold on
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
PlotMuaBoundary(CtrlVar,MUA,'b')  ;
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)

subplot(2,2,2);
QuiverColorGHG(x,y,us-Meas.us,vs-Meas.vs,CtrlVar); axis equal ; title('(us-Meas.us,v-Meas.vs)') ;
hold on ;
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
PlotMuaBoundary(CtrlVar,MUA,'b')  ;
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)

subplot(2,2,3);
[~,~,QuiverPar]=QuiverColorGHG(x,y,Meas.us,Meas.vs,CtrlVar); axis equal ;
title('(Meas.us,Meas.vs)') ;
hold on ;
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
PlotMuaBoundary(CtrlVar,MUA,'b')  ;
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)

subplot(2,2,4);
QuiverPar.QuiverSameVelocityScalingsAsBefore=1;
QuiverColorGHG(x,y,us,vs,QuiverPar); axis equal ; title('(us,vs)') ;
hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
PlotMuaBoundary(CtrlVar,MUA,'b')  ;
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel);
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)

%%  % Difference in speed

SpeedMeas=sqrt(Meas.us.^2+Meas.vs.^2);
SpeedCalc=sqrt(us.^2+vs.^2);

SpeedDiff=100*(SpeedCalc-SpeedMeas)./SpeedMeas;
figure 

PlotMeshScalarVariable(CtrlVar,MUA,SpeedDiff);
hold on 
hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
title('100*(SpeedCalc-SpeedMeas)./SpeedMeas')
axis([min(x) max(x) min(y) max(y)]/CtrlVar.PlotXYscale)

%%

if CtrlVar.Inverse.TestAdjoint.isTrue
    
    dJdp=InvFinalValues.dJdp;
    dJdpTest=InvFinalValues.dJdpTest;
    fprintf('Node/Ele   dJdp          dJdpTest      dJdp-dJdpTest     dJdp/dtdpTest \n')
    iRange=find(~isnan(dJdpTest));
    for ii=1:numel(iRange)
        I=iRange(ii);
        fprintf('%i %15g %15g  %15g  %15g \n',I,dJdp(I),dJdpTest(I),dJdp(I)-dJdpTest(I),dJdp(I)/dJdpTest(I))
    end
    %[dJdp(iRange) dJdpTest(iRange)   dJdp(iRange)-dJdpTest(iRange) dJdp(iRange)./dJdpTest(iRange)]
    fprintf('||dJdpTest-dJdp||/||dJdp||= %g \n ',norm(dJdpTest(iRange)-dJdp(iRange))/norm(dJdp(iRange)))
    
    %%
    
    if ~(isempty(InvFinalValues.dJdC) && isempty(InvFinalValues.dJdCTest))
        
        IFigC=figure('Name','Inversion C','NumberTitle','off');
        
        
        
        subplot(2,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdC) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        title('dJdC Adjoint gradient')
        
        subplot(2,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdCTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        title('dJdC Brute force gradient')
        
        
        subplot(2,2,3) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdC-InvFinalValues.dJdCTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        title('Difference between adjoint and brute force derivatives')
        
        subplot(2,2,4) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdC./InvFinalValues.dJdCTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        title('Ratio between adjoint and brute force derivatives')
        
        IFigC.Position=[948.43 41.571 1246.3 1115.4];
        %%
    end
    
    
    if ~(isempty(InvFinalValues.dJdAGlen) && isempty(InvFinalValues.dJdAGlenTest))
        
        IFigAGlen=figure('Name','Inversion AGlen','NumberTitle','off');
        
        
        
        subplot(2,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlen) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        title('dJdAGlen Adjoint gradient')
        
        subplot(2,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlenTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        title('dJdAGlen Brute force gradient')
        
        
        subplot(2,2,3) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlen-InvFinalValues.dJdAGlenTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        title('Difference between adjoint and brute force derivatives')
        
        subplot(2,2,4) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlen./InvFinalValues.dJdAGlenTest) ;
        hold on
        PlotMuaMesh(CtrlVar,MUA);
        title('Ratio between adjoint and brute force derivatives')
        
        IFigAGlen.Position=[1.5714 41.571 1096 1115.4];
        %%
    end
    
    
    %%
else
    
    if ~isempty(InvFinalValues.dJdAGlen)
        IFigGradientsA=figure('Name','dJdAGlen Gradients','NumberTitle','off');
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdAGlen) ; title('dJdAGlen')
        hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    end
    
    if ~isempty(InvFinalValues.dJdC)
        IFigGradientsC=figure('Name','dJdC Gradients','NumberTitle','off');
        PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.dJdC) ; title('dJdC')
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
            hold on ; PlotMuaMesh(CtrlVar,MUA,[],'w')
            subplot(1,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.C) ; title('Retrieved C')
            hold on ; PlotMuaMesh(CtrlVar,MUA,[],'w')
            tFig1.Units='normalized';
            tFig1.Position=[0.5 0.51 0.5 0.4];
            
            if  ~CtrlVar.CisElementBased
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
        end
    end
    
    if contains(lower(CtrlVar.Inverse.InvertFor),'aglen')
        
        if ~isempty(Priors.TrueAGlen)
            tFig1=figure('Name','True and estimated AGlen','NumberTitle','off');
            subplot(1,2,1) ; PlotMeshScalarVariable(CtrlVar,MUA,Priors.TrueAGlen) ; title('True AGlen')
            hold on ; PlotMuaMesh(CtrlVar,MUA,[],'w')
            subplot(1,2,2) ; PlotMeshScalarVariable(CtrlVar,MUA,InvFinalValues.AGlen) ; title('Retrieved AGlen')
            hold on ; PlotMuaMesh(CtrlVar,MUA,[],'w')
            tFig1.Units='normalized';
            tFig1.Position=[0.5 0.51 0.5 0.4];
            
            if  ~CtrlVar.AGlenisElementBased
                tFig2=figure('Name','Difference between true and estimated','NumberTitle','off');
                %PlotMeshScalarVariable(CtrlVar,MUA,Priors.TrueC-InvFinalValues.C);
                
                subplot(1,3,1)
                tri=PlotNodalVariableAsTriSurface(CtrlVar,MUA,[],InvFinalValues.AGlen);
                xlabel('x') ; ylabel('y') ; title('Inverted AGlen')
                colorbar('south')
                
                subplot(1,3,2)
                PlotNodalVariableAsTriSurface(CtrlVar,MUA,tri,Priors.TrueAGlen);
                xlabel('x') ; ylabel('y') ; title('True AGlen')
                colorbar('south')
                
                subplot(1,3,3)
                PlotNodalVariableAsTriSurface(CtrlVar,MUA,tri,InvFinalValues.AGlen-Priors.TrueAGlen);
                xlabel('x') ; ylabel('y') ; title('True AGlen')
                title('AGlen: True-Estimated')
                colorbar('south')
                
                tFig2.Units='normalized';
                tFig2.Position=[0.1 0.2 0.8 0.5];
            end
        end
    end
    
    
    %%
    if ~isempty(RunInfo.Inverse.J)
        FigObj=figure('Name','Inverse Parameter Optimisation','NumberTitle','off');
        hold off
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
    
    
end
%%
end
