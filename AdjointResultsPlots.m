function AdjointResultsPlots(Experiment,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,alpha,rho,rhow,g,GF,...
                InvStartValues,Priors,Meas,BCsAdjoint,Info,InvFinalValues,xAdjoint,yAdjoint)

%save TestSave
%error('fdas')

%%
if nargin==0
    fprintf('loading AdjointRestart')
    load AdjointRestart
    fprintf(' done.\n')
end

if nargin==1 % if only one input argument is given, it is assumed that it is the name of an input file
    fprintf('loading %s',Experiment)
    load(Experiment)
    fprintf(' done.\n')
end

%%
us=ub+ud; vs=vb+vd;
x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
tri=MUA.connectivity;
GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
xGL=[] ; yGL=[] ; 

figure(3000)

subplot(1,2,1)
[FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,Meas.us,CtrlVar);  title('us Meas on numerical grid') ; colorbar
hold on ; 
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r'); 
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel); 

subplot(1,2,2)
[FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,Meas.vs,CtrlVar);  title('vs Meas on numerical grid') ; colorbar
hold on ; 
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r'); 
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel); 

%%

usError=sqrt(spdiags(Meas.usCov));
vsError=sqrt(spdiags(Meas.vsCov));
wsError=sqrt(spdiags(Meas.wsCov));
figure(3010) 
subplot(2,1,1)
[FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,usError,CtrlVar);  
title('usError on numerical grid') ; colorbar
hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r'); 
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel); 

subplot(2,1,2) 
[FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,vsError,CtrlVar);  
title('vsError on numerical grid') ; colorbar
hold on ; [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r'); 
xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel); 

if strcmp(CtrlVar.AdjointGrad,'A')
    
    figure ; hPatch=PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates/CtrlVar.PlotXYscale,log10(InvFinalValues.AGlen));
    title('log10(InvFinalValues.AGlen)') ; cbar=colorbar; title(cbar, '(a^{-1} kPa^{-3})');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel); 
    
    figure ; hPatch=PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates/CtrlVar.PlotXYscale,log10(InvStartValues.AGlen));
    title('log10(Astart)') ; cbar=colorbar; title(cbar, '(a^{-1} kPa^{-3})');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel); 
    
    figure ; hPatch=PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates/CtrlVar.PlotXYscale,log10(InvFinalValues.AGlen)-log10(InvStartValues.AGlen));
    title('log10(InvFinalValues.AGlen)-log10(InvStartValues.AGlen)') ; cbar=colorbar; title(cbar, '(a^{-1} kPa^{-3})');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel); 
end


if strcmp(CtrlVar.AdjointGrad,'C')
    
    figure ; hPatch=PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates/CtrlVar.PlotXYscale,log10(InvFinalValues.C));
    title('log10(InvFinalValues.C)') ; cbar=colorbar; title(cbar, '(m/a/kPa^m)');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel); 
    
    figure ; hPatch=PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates/CtrlVar.PlotXYscale,log10(InvStartValues.C));
    title('log10(Cstart)') ; cbar=colorbar; title(cbar, '(m/a/kPa^m)');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel); 
    
    figure ; hPatch=PlotElementBasedQuantities(MUA.connectivity,MUA.coordinates/CtrlVar.PlotXYscale,log10(InvFinalValues.C)-log10(InvStartValues.C));
    title('log10(InvFinalValues.C)-log10(Cstart)') ; cbar=colorbar; title(cbar, '(m/a/kPa^m)');
    xlabel(CtrlVar.PlotsXaxisLabel);  ylabel(CtrlVar.PlotsYaxisLabel); 
end

[tbx,tby,tb,beta2] = CalcBasalTraction(CtrlVar,MUA,ub,vb,InvFinalValues.C,InvFinalValues.m,GF) ;
figure
PlotMeshScalarVariable(CtrlVar,MUA,tb) ;
title(' tb ') ; cbar=colorbar; title(cbar, '(kPa)');


CtrlVar.VelPlotIntervalSpacing='log10';
figure(3020)
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
QuiverColorGHG(x,y,Meas.us,Meas.vs,CtrlVar); axis equal ; title('(Meas.us,Meas.vs)') ; 
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
figure   
[It,~]=size(Info.JoptVector);
semilogy(0:It-1,Info.JoptVector(:,1),'-ro') ; hold on
semilogy(0:It-1,Info.JoptVector(:,2),'-b+') ;
semilogy(0:It-1,Info.JoptVector(:,3),'-g*') ;
semilogy(0:It-1,Info.JoptVector(:,5),'-mx') ;
semilogy(0:It-1,Info.JoptVector(:,4),'-cs') ;
semilogy(0:It-1,Info.JoptVector(:,6),'-yd') ;

legend('Cost function','Data misfit','C Reg','C barrier','A Reg','A barrier')
xlabel('Iteration') ;
hold off
%%
end
