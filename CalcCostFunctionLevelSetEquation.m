function [r,UserVar,RunInfo,rForce,rWork,D2]=CalcCostFunctionLevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,gamma,F1,F0,L,Lrhs,l,dLSF,dl)
    
    
    
    narginchk(12,12)
    nargoutchk(1,6)
    
    
    F1.LSF=F1.LSF+gamma*dLSF;
 
    
    [UserVar,R]=LevelSetEquationAssemblyNR2(UserVar,CtrlVar,MUA,F0.LSF,F0.c,F0.ub,F0.vb,F1.LSF,F1.c,F1.ub,F1.vb);
    
    
    if ~isempty(L)
        
        % frhs=-R-L'*l;
        % grhs=Lrhs-L*F1.LSF;
        
        % This is for linear equations which will always be fullfiled 
        frhs=-R/MUA.Area ; 
        grhs=Lrhs-L*F1.LSF;
        
        
    else
        frhs=-R;
        grhs=[];
    end
    
    frhs=frhs/MUA.Area;
    % grhs=grhs;           % This has the units \varphi and there is no intergration over the domain
           
    D2=[frhs;grhs]'*[dLSF;dl]; 
    rWork=full(D2^2);
    
    rForce=full([frhs;grhs]'*[frhs;grhs]); 
    
 
    
    switch CtrlVar.LSFMinimisationQuantity
        case "Force Residuals"
            r=rForce;
        case "Work Residuals"
            r=rWork;
    end
    
    
end

