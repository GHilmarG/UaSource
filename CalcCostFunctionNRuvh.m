function [UserVar,RunInfo,r,rForce,rWork,ruv,rh,D2]=CalcCostFunctionNRuvh(UserVar,RunInfo,CtrlVar,MUA,F1,F0,dub,dvb,dh,dl,L,luvh,cuvh,gamma,fext0)
    
    
    narginchk(15,15)
    nargoutchk(3,8)
    
    
    F1.ub=F1.ub+gamma*dub;
    F1.vb=F1.vb+gamma*dvb;
    F1.h=F1.h+gamma*dh;
    % luvh=luvh+gamma*dl;
    luvh=luvh+dl;
    
    
    
    CtrlVar.uvhMatrixAssembly.ZeroFields=false;
    CtrlVar.uvhMatrixAssembly.Ronly=true;
    
    
    [UserVar,RunInfo,R,~]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1);
    
    
    if ~isempty(L)
        
        frhs=-R-L'*luvh;
        grhs=cuvh-L*[F1.ub;F1.vb;F1.h];
        
    else
        frhs=-R;
        grhs=[];
    end
    
    
    d=[dub;dvb;dh]  ; % Newton step
    
    D2=frhs'*d  ;
    rWork=D2^2 ;
    
    rForce=ResidualCostFunction(CtrlVar,MUA,L,frhs,grhs,fext0,"-uvh-");
    
    
    
    switch CtrlVar.uvhCostFunction
        case "Force Residuals"
            r=rForce;
        case "Work Residuals"
            r=rWork;
    end
    
    if nargout>=6
        Inodes=[]; 
        [ruv,rh]=CalcIncrementsNorm(CtrlVar,MUA,L,Inodes,F1.ub,dub,F1.vb,dvb,F1.h,dh);
    end
    
end

