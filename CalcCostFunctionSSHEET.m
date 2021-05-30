function [r,UserVar,RunInfo,rForce,rWork,D2]=CalcCostFunctionSSHEET(UserVar,RunInfo,CtrlVar,gamma,dh,MUA,AGlen,n,C,m,rho,g,h0,b0,h1,b1,a0,a1,dt,Lh,lh,dlh,fext0,ch)
    
    narginchk(24,24)
    
    %s1=s1+gamma*dh;
    %h=s1-b1;

    h1=h1+gamma*dh; 
    lh=lh+gamma*dlh;
    
 
    R=MatrixAssemblySSHEETtransient2HD(CtrlVar,MUA,AGlen,n,C,m,rho,g,h0,b0,h1,b1,a0,a1,dt);
    
    
    if ~isempty(Lh)
        
        frhs=-R-Lh'*lh;
        grhs=ch-Lh*h1;
        
    else
        frhs=-R;
        grhs=[];
    end
    
    
    rForce=[frhs;grhs]'*[frhs;grhs]./(fext0'*fext0+1000*eps); 
    D2=[frhs;grhs]'*[dh;dlh]  ;
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

