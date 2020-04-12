function [UserVar,RunInfo,r,rRes,rWork,rDisp,D2]=CalcCostFunctionNRuvh(UserVar,RunInfo,CtrlVar,MUA,F1,F0,dub,dvb,dh,dl,L,luvh,cuvh,gamma,fext0)
    
    
    narginchk(15,15)
    nargoutchk(3,7) % later reduce this (TestIng)
    
    % TestIng
    % P=F0.LSF>0 ; 
  
    
    
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
    
    %%
    
    %% TestIng
    
    % frhs=[P;P;P].*frhs;
    % fext0=[P;P;P].*fext0;
    %d=[P;P;P].*d;
    
    d=[dub;dvb;dh]  ; % Newton step
    
    D2=frhs'*d  ;
    rWork=D2^2 ;

   rRes=ResidualCostFunction(CtrlVar,MUA,L,frhs,grhs,fext0,"-uvh-");
        
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

