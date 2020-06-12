function [UserVar,RunInfo,LSF,lambda]=LevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1)
    %%
    %
    %
    %  %  df/dt + u df/dx + v df/dy - div (kappa grad f) = c norm(grad f0)
    %
    
    narginchk(7,7)
    nargoutchk(3,3)
    
    
    persistent LastResetTime
    
    
    if ~CtrlVar.LevelSetMethod
        LSF=F0.LSF;
        lambda=[];
        return
    end
    
    if CtrlVar.CalvingLaw=="-No Ice Shelves-" 
    
        [LSF,UserVar,RunInfo]=ReinitializeLevelSet(UserVar,RunInfo,CtrlVar,MUA,F1.GF.node,CtrlVar.GLthreshold);
        return 
    end
    
    if isempty(LastResetTime)
        LastResetTime=0 ;
    end
 
       
    if CtrlVar.time>( LastResetTime+CtrlVar.LevelSetReinitializeTimeInterval) 
        fprintf("LevelSetEquation: Level Set is re-initialized. \n")
        [F0.LSF,UserVar,RunInfo]=ReinitializeLevelSet(UserVar,RunInfo,CtrlVar,MUA,F0.LSF)  ;
        LastResetTime=CtrlVar.time ; 
    end
    
    
    
    
    switch CtrlVar.LevelSetSolutionMethod
        
        case "Piccard"
            
            [UserVar,RunInfo,LSF,lambda]=LevelSetEquationPiccard(UserVar,RunInfo,CtrlVar,MUA,BCs,F0);
            
        otherwise
            
            % This will actually also do Piccard unless CtrlVar.LevelSetSolutionMethod="Newton-Raphson" ;
            [UserVar,RunInfo,LSF,lambda]=LevelSetEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1);
            
    end
 
    
    
    if CtrlVar.LevelSetInfoLevel>=100 && CtrlVar.doplots
       
        F1.LSF=LSF ; % here needed for plotting
        [fLSF1,fLSF0,fdLSF,fMeshLSF]=LevelSetInfoLevelPlots(CtrlVar,MUA,BCs,F0,F1);
        
    end
    
    
end



