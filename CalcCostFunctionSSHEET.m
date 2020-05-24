function [r,UserVar,RunInfo,rForce,rWork,D2]=CalcCostFunctionSSHEET(UserVar,RunInfo,CtrlVar,gamma,dh,MUA,AGlen,n,rho,g,s0,b0,s1,b1,a0,a1,dt,Lh,lambdah,dlambdah,fext0,ch)
    
    narginchk(22,22)
    
    s1=s1+gamma*dh;
    h=s1-b1;
    R=MatrixAssemblySSHEETtransient2HD(CtrlVar,MUA,AGlen,n,rho,g,s0,b0,s1,b1,a0,a1,dt);
    
    
    if ~isempty(Lh)
        
        frhs=-R-Lh'*(lambdah+gamma*dlambdah);
        grhs=ch-Lh*(h+gamma*dh);
        
    else
        frhs=-R;
        grhs=[];
    end
    
    
    rForce=[frhs;grhs]'*[frhs;grhs]./(fext0'*fext0+1000*eps); 
    D2=[frhs;grhs]'*[dh;dlambdah]  ;
    rWork=D2^2 ;
    
    rForce=full(rForce) ; rWork=full(rWork) ; D2=full(D2); 
    
    
    switch CtrlVar.hMinimisationQuantity
        case "Force Residuals"
            r=rForce;
        case "Work Residuals"
            r=rWork;
        otherwise
            error("CalcCostFunctionNR:UnknownCase")
    end
    
    
end

