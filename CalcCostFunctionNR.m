function [r,UserVar,RunInfo,rForce,rWork,D2,frhs,grhs,Normalisation] = CalcCostFunctionNR(UserVar,RunInfo,CtrlVar,MUA,gamma,F,fext0,L,l,cuv,dub,dvb,dl)
    
    nargoutchk(1,9)
    narginchk(13,13)
    
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
    
    
    % rForce=(frhs'*frhs+grhs'*grhs)/(fext0'*fext0+1000*eps);
    Normalisation=fext0'*fext0+1000*eps;
    rForce=full([frhs;grhs]'*[frhs;grhs]./Normalisation); 
    
    % Newton Decrement
    
    D2=full([frhs;grhs]'*[dub;dvb;dl])  ;
    rWork=D2^2 ;
    
    
    
    switch CtrlVar.uvMinimisationQuantity
        case "Force Residuals"
            r=rForce;
        case "Work Residuals"
            r=rWork;
        otherwise
            error("CalcCostFunctionNR:UnknownCase")
    end
    
    % fprintf('gamma=%f rRes=%g \t rWork=%g \t D2=%g \n',gamma,rRes,rWork,D2)
    
end

