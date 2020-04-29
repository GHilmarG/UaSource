function [UserVar,RunInfo,LSF1,lambda]=LevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,BCsLevelSet,F0,F1)
    %%
    %
    %
    %  %  df/dt + u df/dx + v df/dy - div (kappa grad f) = c norm(grad f0)
    %
    
    persistent iCalls
    
    
    if ~CtrlVar.LevelSetMethod
        LSF1=LSF0;
        lambda=[];
        return
    end
    
    if isempty(iCalls)
        iCalls=0 ;
    end
    iCalls=iCalls+1;
    
    switch CtrlVar.LevelSetSolutionMethod
        
        case "Piccard"
            
            [UserVar,RunInfo,LSF1,lambda]=LevelSetEquationPiccard(UserVar,RunInfo,CtrlVar,MUA,BCsLevelSet,F0); 
            
        otherwise
            
            % This will actually also do Piccard unless CtrlVar.LevelSetSolutionMethod="Newton-Raphson" ;
            [UserVar,RunInfo,LSF1,lambda]=LevelSetEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,BCsLevelSet,F0,F1);
            
            
    end
    
    
    
end

