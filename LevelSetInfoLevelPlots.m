function  [fLSF1,fLSF0,fdLSF,fMeshLSF]=LevelSetInfoLevelPlots(CtrlVar,MUA,BCs,F0,F1)
    
    fLSF1=[] ; fLSF0=[] ; fdLSF=[] ; fMeshLSF=[];
    
    if nargin<5
        F1=[];
    end
    
    
    
    
    
    
    if ~isempty(F1)
        
        fLSF1=FindOrCreateFigure('LSF1');
        hold off ; PlotMeshScalarVariable(CtrlVar,MUA,F1.LSF); title('LSF1')
        hold on ; [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F1,'r');
        colormap(othercolor('BuOr_12',1024));
    end
    
    
    if ~isempty(F0)
        fLSF0=FindOrCreateFigure('LSF0');
        hold off ; PlotMeshScalarVariable(CtrlVar,MUA,F0.LSF); title('LSF0')
        hold on ; [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F0,'r');
        colormap(othercolor('BuOr_12',1024));
    end
    
    if ~isempty(F0) && ~isempty(F1)
        fdLSF=FindOrCreateFigure('dLSF');
        hold off ; PlotMeshScalarVariable(CtrlVar,MUA,F1.LSF-F0.LSF); title('dLSF1')
        hold on ; [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F1,'r');
        colormap(othercolor('BuOr_12',1024));
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
            fMeshLSF.CurrentAxes.Title.String={Temp,"Grounding line in red"};
        end
        if ~isempty(F1.LSF) && CtrlVar.LevelSetMethod
            hold on ; [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F1,'b','LineWidth',2) ;
            Temp=fMeshLSF.CurrentAxes.Title.String;
            fMeshLSF.CurrentAxes.Title.String={Temp,"Calving front in blue"};
        end
        Par.RelativeVelArrowSize=10 ;
        QuiverColorGHG(MUA.coordinates(:,1)/1000,MUA.coordinates(:,2)/1000,F1.ub,F1.vb,Par) ;
        if KeepAxes
            xlim(xL) ; ylim(yL) ;
        end
    end
    
end