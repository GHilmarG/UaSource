function [r,UserVar,RunInfo,rForce,rWork,D2]=CalcCostFunctionLevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,gamma,F1,F0,L,Lrhs,l,dLSF,dl)
    
    
    
    narginchk(12,12)
    nargoutchk(1,6)
    
    
    F1.LSF=F1.LSF+gamma*dLSF;
    l=l+gamma*dl;
    
%     % Testing
%     F0.LSFqx=zeros(MUA.Nnodes,1) ;
%     F0.LSFqy=zeros(MUA.Nnodes,1) ;
%     F1.LSFqx=zeros(MUA.Nnodes,1) ;
%     F1.LSFqy=zeros(MUA.Nnodes,1) ;
%     % 
    
    [UserVar,R]=LevelSetEquationAssemblyNR2(UserVar,CtrlVar,MUA,F0.LSF,F0.c,F0.ub,F0.vb,F1.LSF,F1.c,F1.ub,F1.vb,F0.LSFqx,F0.LSFqy,F1.LSFqx,F1.LSFqy);
    
    
    if ~isempty(L)
     
         frhs=-R-L'*l;        % Units: Area  [\varphi] / time 
         grhs=Lrhs-L*F1.LSF;  % Units: [\varphi]=distance, but this should always be zero anyhow if initial point is feasable
        
        
    else
        frhs=-R;
        grhs=[];
    end
    
    frhs=frhs/MUA.Area;  % It's OK to do this only here, because I scale frhs and grhs equally.
    grhs=grhs/MUA.Area;
    
           
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

