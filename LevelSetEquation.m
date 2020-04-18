function [UserVar,phi1,lambda]=LevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,BCsLevelSet,F,phi0)
    %%
    %
    %
    %  %  df/dt + u df/dx + v df/dy - div (kappa grad f) = c norm(grad f0)
    %
    
    persistent iCalls
    
    
    if ~CtrlVar.LevelSetMethod
        phi1=phi0; 
        lambda=[]; 
        return
    end
    
    
    
    if isempty(iCalls)
        iCalls=0 ;
    end
    iCalls=iCalls+1;
    
    
    % Define calving rate
    % F.c=zeros(MUA.Nnodes,1)-100e3 ;
    
    
    if ~isfield(CtrlVar,'LevelSetResetInterval') || isempty(CtrlVar.LevelSetResetInterval)
        CtrlVar.LevelSetResetInterval=10000;
    end
    
    MLC=BCs2MLC(CtrlVar,MUA,BCsLevelSet);
    L=MLC.hL ; Lrhs=MLC.hRhs ;
    
    
    kappa=0 ;  % for the time being
    [UserVar,kv,rh]=LevelSetEquationAssembly(UserVar,CtrlVar,MUA,phi0,F.c,F.ub,F.vb,kappa);
    [phi1,lambda]=solveKApe(kv,L,rh,Lrhs,[],[],CtrlVar);
    phi1=full(phi1);
    
    
    if isempty(phi0) ||  mod(iCalls,CtrlVar.LevelSetResetInterval)==0
        
        % >>>>>  Initialize
        CtrlVar.PlotGLs=0 ;
        GF.node=phi1;
        [xc,yc]=PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'k') ;
        
        x=MUA.coordinates(:,1);
        y=MUA.coordinates(:,2);
        xMax=max(x) ; xMin=min(x) ;  xL=xMax-xMin;
        yMax=max(y) ; yMin=min(y) ;  yL=yMax-yMin;
        CalvingFrontClosure=[xc(end) yMax+yL ;   xMin-xL yMax+yL ; xMin-xL yMin-yL ;  xc(end) yMin-yL ] ;
        
        CalvingFront=[xc(:) yc(:) ; CalvingFrontClosure ] ;
        Npoints=1000 ; CalvingFront = interparc(Npoints,CalvingFront(:,1),CalvingFront(:,2),'linear'); % add some points
        DistSigned=SignedDistance(MUA.coordinates,CalvingFront);
        phi1=DistSigned ;
        % <<<<<<
    end
    
    
    
    
    
    %     %% Plotting
    %     figLSF=FindOrCreateFigure("Level Set Field");
    %     PlotMeshScalarVariable(CtrlVar,MUA,phi1) ; title("Level Set")
    %     hold on
    %     F.LSF=phi1 ;
    %     [xc,yc]=PlotCalvingFronts(CtrlVar,MUA,F,'r');
    
end

