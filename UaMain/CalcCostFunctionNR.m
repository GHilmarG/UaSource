




function [r,UserVar,RunInfo,rForce,rWork,D2,frhs,grhs,Normalisation] = CalcCostFunctionNR(UserVar,RunInfo,CtrlVar,MUA,gamma,F,fext0,L,l,cuv,dub,dvb,dl)
    
    nargoutchk(1,9)
    narginchk(13,13)
    
    if isnan(gamma) ; error(' gamma is nan ') ; end
    if ~isreal(gamma) ; error(' gamma is not real ') ; end
    
    F.ub=F.ub+gamma*dub;
    F.vb=F.vb+gamma*dvb;
    l.ubvb=l.ubvb+gamma*dl;
    
   
    CtrlVar.uvMatrixAssembly.Ronly=true;  
    % Ruv=KRTFgeneralBCs(CtrlVar,MUA,F);
    [RunInfo,Ruv]=uvMatrixAssembly(RunInfo,CtrlVar,MUA,F);
    
    if ~isempty(L)
        frhs=-Ruv-L'*l.ubvb;
        grhs=cuv-L*[F.ub;F.vb];
        
    else
        frhs=-Ruv;
        grhs=[];
        dl=[];
    end
    
    
    % (15 Oct 2026) The residual measures are calculated in uvResidualMeasures.m (same code as previously here). The force residual
    % is always the merit function (the option CtrlVar.uvMinimisationQuantity="Work Residuals" has been removed).
    [rForce,rWork,D2,Normalisation]=uvResidualMeasures(frhs,grhs,fext0,dub,dvb,dl) ;
    r=rForce ;
    
    % fprintf('gamma=%f rRes=%g \t rWork=%g \t D2=%g \n',gamma,rRes,rWork,D2)
    
end

