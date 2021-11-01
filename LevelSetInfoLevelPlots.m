function  [fLSF1,fLSF0,fdLSF,fMeshLSF]=LevelSetInfoLevelPlots(CtrlVar,MUA,BCs,F0,F1)
    
    fLSF1=[] ; fLSF0=[] ; fdLSF=[] ; fMeshLSF=[];
    
    if nargin<5
        F1=[];
    end
    
    
    figBC=FindOrCreateFigure("BCs and LSF") ; 
    hold off
    lgd=PlotBoundaryConditions(CtrlVar,MUA,BCs) ;
    hold on ;
    CtrlVar.LineUpGLs=true ; 
    [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F1,LineWidth=2,color="r") ;
    plot(xc/CtrlVar.PlotXYscale,yc/CtrlVar.PlotXYscale,'*r',DisplayName="C") ;
    lgd.String{1}='BCs for $\varphi$';
    lgd.String{2}='$\varphi$';
    lgd.String{2}='$(x_c,y_c)$';
    xlabel("$x\,(\mathrm{km})$",interpreter="latex")
    ylabel("$y\,(\mathrm{km})$",interpreter="latex")
    
    FindOrCreateFigure("F1.LSF") ; hold off
    PlotMeshScalarVariable(CtrlVar,MUA,F1.LSF) ; 
    hold on ; 
    PlotCalvingFronts(CtrlVar,MUA,F1,'r') ;
    title(sprintf("%s F1.LSF at time=%5.2f ",CtrlVar.LevelSetPhase,CtrlVar.time))
    
    FindOrCreateFigure("F1.LSF surface") ; 
    hold off
    PatchObject=PlotMeshScalarVariableAsSurface(CtrlVar,MUA,F1.LSF/CtrlVar.PlotXYscale,500) ; 
    colormap(othercolor('BuOr_12',1024));
    zlabel("$\varphi_1$",interpreter="latex")
    title(sprintf("%s F1.LSF at time=%5.2f ",CtrlVar.LevelSetPhase,CtrlVar.time))
    
    FigC=FindOrCreateFigure("Calving rate");
    PlotMeshScalarVariable(CtrlVar,MUA,F1.c) ;
    xlabel("$x\,(\mathrm{km})$",interpreter="latex")
    ylabel("$y\,(\mathrm{km})$",interpreter="latex")
    title(sprintf("Calving rate c at time=%5.2f ",CtrlVar.time))
    
    speed=sqrt(F1.ub.^2+F1.vb.^2);
    FigC=FindOrCreateFigure("speed");
    PlotMeshScalarVariable(CtrlVar,MUA,speed) ;
    xlabel("$x\,(\mathrm{km})$",interpreter="latex")
    ylabel("$y\,(\mathrm{km})$",interpreter="latex")
    title(sprintf("speed at time=%5.2f ",CtrlVar.time))
    
    
    FigC=FindOrCreateFigure("speedc");
    PlotMeshScalarVariable(CtrlVar,MUA,speed+F1.c) ;
    xlabel("$x\,(\mathrm{km})$",interpreter="latex")
    ylabel("$y\,(\mathrm{km})$",interpreter="latex")
    title(sprintf("speed+c at time=%5.2f ",CtrlVar.time))
    
    
    if ~isempty(F1)
        
        fLSF1=FindOrCreateFigure('LSF1');
        hold off ; 
        [~,ColorbarHandle]=PlotMeshScalarVariable(CtrlVar,MUA,F1.LSF/CtrlVar.PlotXYscale); title('LSF1')
        hold on ; 
        hold on ; [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F1,'k',LineWidth=2);
        colormap(othercolor('BuOr_12',1024));
        title(sprintf("%s F1.LSF at time=%5.2f ",CtrlVar.LevelSetPhase,CtrlVar.time))
        xlabel("$x\,(\mathrm{km})$",interpreter="latex")
        ylabel("$y\,(\mathrm{km})$",interpreter="latex")
        title(ColorbarHandle,"$\varphi_1\,\mathrm{(km)}$",interpreter="latex")
        axis tight
    end
    
    
    if ~isempty(F0)
        fLSF0=FindOrCreateFigure('LSF0');
        hold off ; 
        [~,ColorbarHandle]=PlotMeshScalarVariable(CtrlVar,MUA,F0.LSF/CtrlVar.PlotXYscale); 
        hold on ; [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F0,'k',LineWidth=2);
        colormap(othercolor('BuOr_12',1024));
        title(sprintf("%s F0.LSF at time=%5.2f ",CtrlVar.LevelSetPhase,CtrlVar.time))
        title(ColorbarHandle,"$\varphi_0\, \mathrm{(km)}$",interpreter="latex")
        xlabel("$x\,(\mathrm{km})$",interpreter="latex")
        ylabel("$y\,(\mathrm{km})$",interpreter="latex")
        axis tight
    end
    
    if ~isempty(F0) && ~isempty(F1)
        fdLSF=FindOrCreateFigure('dLSF');
        hold off ; 
        [~,ColorbarHandle]=PlotMeshScalarVariable(CtrlVar,MUA,(F1.LSF-F0.LSF)/CtrlVar.PlotXYscale); title('dLSF1')
        hold on ; [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F1,'k',LineWidth=2);
        colormap(othercolor('BuOr_12',1024));
        title(sprintf("%s F1.LSF-F0.LSF at time=%5.2f ",CtrlVar.LevelSetPhase,CtrlVar.time))
        xlabel("$x\,(\mathrm{km})$",interpreter="latex")
        ylabel("$y\,(\mathrm{km})$",interpreter="latex")
        title(ColorbarHandle,"$\varphi_1-\varphi_0$",interpreter="latex")
        axis tight
    end
    
    if ~isempty(F1)
        fMeshLSF=FindOrCreateFigure("Mesh and LSF");
        if ~isempty(fMeshLSF.CurrentAxes)
            KeepAxes=true;
            xL=fMeshLSF.CurrentAxes.XLim ; yL=fMeshLSF.CurrentAxes.YLim ;
        else
            KeepAxes=false;
        end
        clf(fMeshLSF) ;
        
        PlotMuaMesh(CtrlVar,MUA); hold on
        [xGL,yGL]=PlotGroundingLines(CtrlVar,MUA,F1.GF,[],[],[],'r','LineWidth',2);
        if ~isempty(xGL)
            Temp=fMeshLSF.CurrentAxes.Title.String;
            fMeshLSF.CurrentAxes.Title.String=[Temp(:)',{"Grounding line in red"}];
        end
        if ~isempty(F1.LSF) && CtrlVar.LevelSetMethod
            hold on ; [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F1,'b','LineWidth',2) ;
            Temp=fMeshLSF.CurrentAxes.Title.String;
            fMeshLSF.CurrentAxes.Title.String=[Temp(:)',{"Calving front in blue"}];
        end
        Par.RelativeVelArrowSize=10 ;
        QuiverColorGHG(MUA.coordinates(:,1)/1000,MUA.coordinates(:,2)/1000,F1.ub,F1.vb,Par) ;
        if KeepAxes
            xlim(xL) ; ylim(yL) ;
        end
        title(sprintf("%s time=%5.2f ",CtrlVar.LevelSetPhase,CtrlVar.time))
    end
    
    
    fMeshLSF=FindOrCreateFigure("F0: norm of slope ");
    [dfdx,dfdy]=calcFEderivativesMUA(F0.LSF,MUA,CtrlVar);
    N0=sqrt(dfdx.*dfdx+dfdy.*dfdy) ;
    
    [FigHandle,ColorbarHandle]=PlotMeshScalarVariable(CtrlVar,MUA,N0) ;
    
    
    fMeshLSF2=FindOrCreateFigure("F1: norm of slope ");
    [dfdx,dfdy]=calcFEderivativesMUA(F1.LSF,MUA,CtrlVar);
    N1=sqrt(dfdx.*dfdx+dfdy.*dfdy) ;
    
    [FigHandle,ColorbarHandle]=PlotMeshScalarVariable(CtrlVar,MUA,N1) ;
    
    
    [N0,N1]=ProjectFintOntoNodes(MUA,N0,N1);
    
    fMeshLSF2=FindOrCreateFigure("F0 node: norm of slope ");
    [FigHandle,ColorbarHandle]=PlotMeshScalarVariable(CtrlVar,MUA,N0 ) ;
    colormap(othercolor('BuOr_12',1024));
    hold on
    hold on ; [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F0,'k',LineWidth=2);
    title(ColorbarHandle,"$\| \nabla \varphi_0 \|$",interpreter="latex")
    xlabel("$x\,(\mathrm{km})$",interpreter="latex")
    ylabel("$y\,(\mathrm{km})$",interpreter="latex")
    axis tight 
    
    fMeshLSF2=FindOrCreateFigure("F1 node: norm of slope ");
    [FigHandle,ColorbarHandle]=PlotMeshScalarVariable(CtrlVar,MUA,N1) ;
    
    colormap(othercolor('BuOr_12',1024));
    hold on
    hold on ; [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F1,'k',LineWidth=2);
    title(ColorbarHandle,"$\| \nabla \varphi_1 \|$",interpreter="latex")
    xlabel("$x\,(\mathrm{km})$",interpreter="latex")
    ylabel("$y\,(\mathrm{km})$",interpreter="latex")
    axis tight
    
end