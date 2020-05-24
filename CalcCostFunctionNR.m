function [UserVar,r,rForce,rWork,D2] = CalcCostFunctionNR(UserVar,CtrlVar,MUA,gamma,F,fext0,L,l,cuv,dub,dvb,dl)
    
    nargoutchk(2,5)
    narginchk(12,12)
    
    if isnan(gamma) ; error(' gamma is nan ') ; end
    if ~isreal(gamma) ; error(' gamma is not real ') ; end
    
    F.ub=F.ub+gamma*dub;
    F.vb=F.vb+gamma*dvb;
    l.ubvb=l.ubvb+gamma*dl;
    
    Ruv=KRTFgeneralBCs(CtrlVar,MUA,F);
    
    if ~isempty(L)
        frhs=-Ruv-L'*l.ubvb;
        grhs=cuv-L*[F.ub;F.vb];
        
    else
        frhs=-Ruv;
        grhs=[];
        dl=[];
    end
    
        
    rForce=(frhs'*frhs+grhs'*grhs)/(fext0'*fext0+1000*eps); 
    
    % Newton Decrement
    
    D2=[frhs;grhs]'*[dub;dvb;dl]  ;
    rWork=D2^2 ;

    
    
    switch CtrlVar.uvMinimisationQuantity
        case "Force Residuals"
            r=full(rForce);
        case "Work Residuals"
            r=full(rWork);
    end
    
    % fprintf('gamma=%f rRes=%g \t rWork=%g \t D2=%g \n',gamma,rRes,rWork,D2) 
    
end

