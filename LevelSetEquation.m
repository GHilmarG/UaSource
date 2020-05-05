function [UserVar,RunInfo,LSF,lambda]=LevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,BCsLevelSet,F0,F1)
    %%
    %
    %
    %  %  df/dt + u df/dx + v df/dy - div (kappa grad f) = c norm(grad f0)
    %
    
    narginchk(7,7)
    nargoutchk(3,3)
    
    
    persistent iCalls
    
    
    if ~CtrlVar.LevelSetMethod
        LSF=F0.LSF;
        lambda=[];
        return
    end
    
    if isempty(iCalls)
        iCalls=0 ;
    end
    iCalls=iCalls+1;
    
    
    
    
    switch CtrlVar.LevelSetSolutionMethod
        
        case "Piccard"
            
            [UserVar,RunInfo,LSF,lambda]=LevelSetEquationPiccard(UserVar,RunInfo,CtrlVar,MUA,BCsLevelSet,F0);
            
        otherwise
            
            % This will actually also do Piccard unless CtrlVar.LevelSetSolutionMethod="Newton-Raphson" ;
            [UserVar,RunInfo,LSF,lambda]=LevelSetEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCsLevelSet,F0,F1);
            
    end
    
    if mod(iCalls,CtrlVar.LevelSetMethodReinitializeInterval)==0
        fprintf("LevelSetEquation: Level Set is re-initialized. \n")
        [LSF,UserVar,RunInfo]=ReinitializeLevelSet(UserVar,RunInfo,CtrlVar,MUA,LSF)  ;
    end
    
    
end

