function [UserVar,RunInfo,r,rRes,rWork,rDisp,D2]=CalcCostFunctionNRuvh(UserVar,RunInfo,CtrlVar,MUA,F1,F0,dub,dvb,dh,dl,L,luvh,cuvh,gamma,fext0)
    
  
    narginchk(15,15)
    nargoutchk(3,7) % later reduce this (TestIng)
    
    
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
    
    rRes=ResidualCostFunction(CtrlVar,MUA,L,frhs,grhs,fext0,"-uvh-");
    
    D2=frhs'*[dub;dvb;dh]  ;
    rWork=D2^2 ;
    
    rDisp=NaN ;
    
    
    
    switch CtrlVar.uvhCostFunction
        case "Force Residuals"
            r=rRes;
        case "Work Residuals"
            r=rWork;
        case "Increments"
            r=rDisp ;
    end
    
    
    
end

