function [r,UserVar,RunInfo,rForce,rWork,D2]=CalcCostFunctionLevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,gamma,F1,F0,L,l,dLSF,dl)
    
    
    
    narginchk(11,11)
    nargoutchk(1,6)
    
    
    F1.LSF=F1.LSF+gamma*dLSF;
    l=l+gamma*dl;

    
    [UserVar,R]=LevelSetEquationAssemblyNR(UserVar,CtrlVar,MUA,F0.LSF,F0.c,F0.ub,F0.vb,F1.LSF,F1.c,F1.ub,F1.vb);
    
    
    if ~isempty(L)
        
        frhs=-R-L'*l;  
        grhs=cuvh-L*F1.LSF;
        
    else
        frhs=-R;
        grhs=[];
    end
    
    frhs=frhs/MUA.Area;
    grhs=grhs/MUA.Area; 
    
    D2=[frhs;grhs]'*[dLSF;dl]; 
    rWork=D2^2;
    
    rForce=full([frhs;grhs]'*[frhs;grhs]); 
    
    switch CtrlVar.LSFMinimisationQuantity
        case "Force Residuals"
            r=rForce;
        case "Work Residuals"
            r=rWork;
    end
    
    
end

