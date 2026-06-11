
function PFFPlots(UserVar,CtrlVar,MUA,F,BCs,BCsphi,phi,Psi,e,PlotTitle,options) 




arguments
    UserVar struct
    CtrlVar struct
    MUA     struct
    F       {mustBeA(F,{'struct','UaFields','numeric'})}
    BCs     {mustBeA(BCs,{'struct','BoundaryConditions','numeric'})}
    BCsphi  {mustBeA(BCsphi,{'struct','BoundaryConditions','numeric'})}
    phi     (:,1) {mustBeNumeric}
    Psi     (:,1) {mustBeNumeric}
    e       (:,1) {mustBeNumeric}
    PlotTitle string=""
    options.CreateVideos logical=false
    options.CloseVideos logical=false

end



persistent phiVideo MeshVideo uvVideo eVideo phiLast

narginchk(10,10)


if ~isfield(CtrlVar.PhaseFieldFracture,"Video")
    CtrlVar.PhaseFieldFracture.Video=false;
end

if ~isfield(UserVar,"VideoFileName")
    UserVar.VideoFileName="";
end

if ~isfield(UserVar,"Experiment")
    UserVar.Experiment="";
end

[xphi,yphi]=CalcMuaFieldsContourLine(CtrlVar,MUA,phi,0.8) ;


figBCs=FindOrCreateFigure("BCs") ; clf(figBCs) ;
PlotBoundaryConditions(CtrlVar,MUA,BCs);
hold on ; plot(xphi/CtrlVar.PlotXYscale,yphi/CtrlVar.PlotXYscale,Color="r",LineWidth=2)
title(PlotTitle)
xlabel("$x$ (km)",Interpreter="latex") ; ylabel("$y$ (km)",Interpreter="latex") ;

FindOrCreateFigure("BCs Phi") ; PlotBoundaryConditions(CtrlVar,MUA,BCsphi);
hold on ; plot(xphi/CtrlVar.PlotXYscale,yphi/CtrlVar.PlotXYscale,Color="r",LineWidth=2)
xlabel("$x$ (km)",Interpreter="latex") ; ylabel("$y$ (km)",Interpreter="latex") ;

if CtrlVar.PhaseFieldFracture.Formulation=="-elastic-" 
FigTitle="Effective Shear Modulus (G)";
else
    FigTitle="Effective rate factor (A)";
end

UaPlots(CtrlVar,MUA,F,F.AGlen,FigureTitle=FigTitle) ; 


set(gca,'ColorScale','log')
hold on ;  plot(xphi/CtrlVar.PlotXYscale,yphi/CtrlVar.PlotXYscale,Color="r",LineWidth=2)
title(sprintf("$A$ ")+PlotTitle,Interpreter="latex")
xlabel("$x$ (km)",Interpreter="latex") ; ylabel("$y$ (km)",Interpreter="latex") ;


%% uv 
uvVideoFile="uv-"+UserVar.Experiment+UserVar.VideoFileName;
if options.CreateVideos
    if isempty(uvVideo)
        uvVideo=VideoWriter(uvVideoFile,"MPEG-4") ;
        uvVideo.FrameRate=1;
        open(uvVideo)
    end
end

fvel=FindOrCreateFigure("uv") ; clf(fvel)
QuiverColorGHG(F.x,F.y,F.ub,F.vb,CtrlVar) ;
hold on ; 
plot(xphi/CtrlVar.PlotXYscale,yphi/CtrlVar.PlotXYscale,Color="r",LineWidth=2)

Tphi=sprintf("$(u,v)$ with $l=$%g m, $G_c$=%g",CtrlVar.PhaseFieldFracture.l,CtrlVar.PhaseFieldFracture.Gc)   ; 
title(Tphi,Interpreter="latex")
subtitle(PlotTitle,Interpreter="latex");
hold on ; PlotMuaBoundary(CtrlVar,MUA,"b")
axis tight 
%vel.Position=[900 70 1200 1200]; 

if options.CreateVideos
    CurFig=gcf; CurFig.Position=[900 70 1200 1200];
    if options.CloseVideos 
        close(uvVideo)
    else
        frame=getframe(gcf);
        writeVideo(uvVideo,frame) ;
    end
end
%% phi

phiVideoFile="phi-"+UserVar.Experiment+UserVar.VideoFileName;


if options.CreateVideos
    if isempty(phiVideo)
        phiVideo=VideoWriter(phiVideoFile,"MPEG-4")  ;
        phiVideo.FrameRate=1;
        open(phiVideo)
    end
end

figphi=FindOrCreateFigure("phi")  ; clf(figphi) ;
cbar=UaPlots(CtrlVar,MUA,F,phi,GroundingLineColor="k") ;
%CM=cmocean('-ice',15) ; colormap(CM);
CL=clim;
% if min(CL) < 0.1 && max(CL)>0.6
%     CM=cmocean('balanced',25,'pivot',0.1) ; colormap(CM);
% else
%     CM=cmocean('balanced',25) ; colormap(CM);
% end
title(cbar,"$\phi$",interpreter="latex")
xlabel("$x$ (km)",Interpreter="latex") ; ylabel("$y$ (km)",Interpreter="latex") ;
clim([0 1])
% pivot=0.5 ;
% if min(phi) < pivot && max(phi) > pivot
% 
%     CM=cmocean('balanced',25,'pivot',pivot) ; colormap(CM);
% 
% end

hold on ;  plot(xphi/CtrlVar.PlotXYscale,yphi/CtrlVar.PlotXYscale,Color="r",LineWidth=2)

Tphi=sprintf("Phase field, $\\phi$, with $l=$%g m, $G_c$=%g",CtrlVar.PhaseFieldFracture.l,CtrlVar.PhaseFieldFracture.Gc)   ; 
title(Tphi,Interpreter="latex")
subtitle(PlotTitle,Interpreter="latex");
xlabel("$x$ (km)",Interpreter="latex") ; ylabel("$y$ (km)",Interpreter="latex")
CM=cmocean('-ice',15) ; colormap(CM);
if options.CreateVideos

    CurFig=gcf; CurFig.Position=[25 70 1200 1200]; axis tight
    if options.CloseVideos 
        close(phiVideo)
    else
        frame=getframe(gcf);
        writeVideo(phiVideo,frame) ;
    end
end
%% Mesh


MeshVideoFile="Mesh-"+UserVar.Experiment+UserVar.VideoFileName;

if options.CreateVideos
    if isempty(MeshVideo)
        MeshVideo=VideoWriter(MeshVideoFile,"MPEG-4") ;
        MeshVideo.FrameRate=1;
        open(MeshVideo)
    end
end


figMesh=FindOrCreateFigure("Mesh")  ; clf(figMesh) ;
PlotMuaMesh(CtrlVar,MUA) ;
axis tight
hold on ;  plot(xphi/CtrlVar.PlotXYscale,yphi/CtrlVar.PlotXYscale,Color="r",LineWidth=2)
title(sprintf("Mesh ")+PlotTitle,Interpreter="latex")
xlabel("$x$ (km)",Interpreter="latex") ; ylabel("$y$ (km)",Interpreter="latex")


if options.CreateVideos

    CurFig=gcf; CurFig.Position=[25 70 1200 1200]; axis tight


    if options.CloseVideos

        close(MeshVideo)

    else

        frame=getframe(gcf);
        writeVideo(MeshVideo,frame) ;

    end
end

%%



% figphiy=FindOrCreateFigure("Phi(y)") ; clf(figphiy) ; 
% Ind=F.x>50e3 & F.x <60e3 ;   
% plot(F.y(Ind)/CtrlVar.PlotXYscale,phi(Ind),'.r') ;
% 
% 
% 
% cbar=UaPlots(CtrlVar,MUA,F,F.rho,FigureTitle="rho effective") ;
% title("$\rho$ effective "+PlotTitle,interpreter="latex")
% title(cbar,"$\rho$",interpreter="latex")
% hold on ;  plot(xphi/CtrlVar.PlotXYscale,yphi/CtrlVar.PlotXYscale,Color="r",LineWidth=2)
% xlabel("$x$ (km)",Interpreter="latex") ; ylabel("$y$ (km)",Interpreter="latex") ;
% 
% cbar=UaPlots(CtrlVar,MUA,F,F.h,FigureTitle="ice thickness (h)") ;
% title("ice thickness ($h$) "+PlotTitle,interpreter="latex")
% title(cbar,"$h$",interpreter="latex")
% hold on ;  plot(xphi/CtrlVar.PlotXYscale,yphi/CtrlVar.PlotXYscale,Color="r",LineWidth=2)
% xlabel("$x$ (km)",Interpreter="latex") ; ylabel("$y$ (km)",Interpreter="latex") ;

%%
eVideoFile="e-"+UserVar.Experiment+UserVar.VideoFileName;
if options.CreateVideos
    if isempty(eVideo)
        eVideo=VideoWriter(eVideoFile,"MPEG-4") ; 
        eVideo.FrameRate=1;
        open(eVideo)
    end
end
fige=FindOrCreateFigure("e") ;  clf(fige);
cbar=UaPlots(CtrlVar,MUA,F,e);
hold on ;
plot(xphi/CtrlVar.PlotXYscale,yphi/CtrlVar.PlotXYscale,Color="r",LineWidth=2)
title(sprintf("Effective strain rate, $\\dot{e}$  (1/yr)")+PlotTitle,Interpreter="latex")
title(cbar,sprintf("$\\dot{e}$  "),Interpreter="latex")
xlabel("$x$ (km)",Interpreter="latex") ; ylabel("$y$ (km)",Interpreter="latex") ;

if CtrlVar.PhaseFieldFracture.Video
    CurFig=gcf; CurFig.Position=[25 70 1200 1200];
    if CtrlVar.PhaseFieldFracture.iphiUpdate==CtrlVar.PhaseFieldFracture.MaxUpdates
        close(eVideo)
    else
        frame=getframe(gcf);
        writeVideo(eVideo,frame) ;
    end
end


cbar=UaPlots(CtrlVar,MUA,F,Psi,FigureTitle="Psi") ;  set(gca,'ColorScale','log')
title(cbar,"$\Psi$",interpreter="latex")
hold on ;  plot(xphi/CtrlVar.PlotXYscale,yphi/CtrlVar.PlotXYscale,Color="r",LineWidth=2)
title(sprintf("Strain energy density, $\\Psi$  ")+PlotTitle,Interpreter="latex")
xlabel("$x$ (km)",Interpreter="latex") ; ylabel("$y$ (km)",Interpreter="latex") ;


cbar=UaPlots(CtrlVar,MUA,F,[F.ub F.vb],FigureTitle="elastic displacements") ;  
title(cbar,"displacements",interpreter="latex")
hold on ;  plot(xphi/CtrlVar.PlotXYscale,yphi/CtrlVar.PlotXYscale,Color="r",LineWidth=2)
title(sprintf("Elastic displacement vectors")+PlotTitle,Interpreter="latex")
xlabel("$x$ (km)",Interpreter="latex") ; ylabel("$y$ (km)",Interpreter="latex") ;


cbar=UaPlots(CtrlVar,MUA,F,F.h,FigureTitle="h",PlotUnderMesh=true) ;  
title(cbar,"$h$ (m)",interpreter="latex")
hold on ;  plot(xphi/CtrlVar.PlotXYscale,yphi/CtrlVar.PlotXYscale,Color="r",LineWidth=2)
title(sprintf("Ice thickness, $h$  ")+PlotTitle,Interpreter="latex")
xlabel("$x$ (km)",Interpreter="latex") ; ylabel("$y$ (km)",Interpreter="latex") ;

% only works if same mesh
% if ~isempty(phiLast)
%     cbar=UaPlots(CtrlVar,MUA,F,phi-phiLast,FigureTitle="dphi",PlotUnderMesh=true) ;
%     title(cbar,"$\Delta \phi$ (m)",interpreter="latex")
%     hold on ;  plot(xphi/CtrlVar.PlotXYscale,yphi/CtrlVar.PlotXYscale,Color="r",LineWidth=2)
%     title(sprintf("Change in phase field, $\\Delta \\phi$  ")+PlotTitle,Interpreter="latex")
%     xlabel("$x$ (km)",Interpreter="latex") ; ylabel("$y$ (km)",Interpreter="latex") ;
% end


phiLast=phi; 


return

% deviatoric stresses : This takes some time
[X,Y]=ndgrid(linspace(min(F.x),max(F.x),80),linspace(min(F.y),max(F.y),80));
I=nearestNeighbor(MUA.TR,[X(:) Y(:)]);  % find nodes within computational grid closest to the regularly scape X and Y grid points.
fstress=FindOrCreateFigure("dev stresses") ; clf(fstress)
scale=5e-3;

iphi=phi>0.5 ;

txx(iphi)=nan;
txy(iphi)=nan;
tyy(iphi)=nan;

PlotTensor(F.x(I)/CtrlVar.PlotXYscale,F.y(I)/CtrlVar.PlotXYscale,txx(I),txy(I),tyy(I),scale);
hold on ; 
plot(xphi/CtrlVar.PlotXYscale,yphi/CtrlVar.PlotXYscale,Color="r",LineWidth=2)
PlotMuaBoundary(CtrlVar,MUA,'k') ; axis equal





end