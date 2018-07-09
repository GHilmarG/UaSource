function Hessian=fminconHessianFcn(x,lambda,MUA,AGlen,C,CtrlVar)

persistent PHessian

if isempty(PHessian)
    
    if (strcmpi(CtrlVar.AdjointGrad,'C') &&  CtrlVar.CisElementBased) || (strcmpi(CtrlVar.AdjointGrad,'A') &&  CtrlVar.AGlenisElementBased)
        
        PHessian=sparse(1:MUA.Nele,1:MUA.Nele,1);
    else
        PHessian=MassMatrix2D1dof(MUA);
        
    end
end

Hessian=PHessian;

end