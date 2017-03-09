function [J,dJdp,Hessian,JGHouts,F]=JGH(p,UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo)

% Calculates objective function, gradient (accurate), Hessian (guessed)

persistent ubP vbP

if CtrlVar.Inverse.ResetPersistentVariables
    ubP=[];
    vbP=[];
end

if ~isempty(ubP)
    F.ub=ubP;
    F.vb=vbP;
end

%
% if CtrlVar.OnlyGetPersistenValues
%     J=[] ; Gradient=[] ;Hessian=[] ; Idata=[] ; IReg=[] ; IBarrier=[] ; l=[];
%     return
%     
% end



switch nargout
    
    case 1
        
        CtrlVar.Inverse.CalcGradR=0;
        CtrlVar.Inverse.CalcGradI=0;
        CtrlVar.Inverse.CalcHessR=0;
        CtrlVar.Inverse.CalcHessI=0;
        
    case 2
        
        CtrlVar.Inverse.CalcGradR=1;
        CtrlVar.Inverse.CalcGradI=1;
        CtrlVar.Inverse.CalcHessR=0;
        CtrlVar.Inverse.CalcHessI=0;
        
    case {3,4,5}
        
        CtrlVar.Inverse.CalcGradR=1;
        CtrlVar.Inverse.CalcGradI=1;
        CtrlVar.Inverse.CalcHessR=1;
        CtrlVar.Inverse.CalcHessI=1;
        
end


if contains(lower(CtrlVar.Inverse.InvertFor),'aglen') && contains(lower(CtrlVar.Inverse.InvertFor),'c')  % AC
    
    NA=numel(InvStartValues.AGlen);
    
    if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
        F.AGlen=10.^p(1:NA);
    else
        F.AGlen=p(1:NA);
    end
    
    if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
        F.C=10.^p(NA+1:end);
    else
        F.C=p(NA+1:end);
    end
    
elseif contains(lower(CtrlVar.Inverse.InvertFor),'aglen')   % A
    
    
    if contains(lower(CtrlVar.Inverse.InvertFor),'logaglen')
        F.AGlen=10.^p;
    else
        F.AGlen=p;
    end
    
    
elseif contains(lower(CtrlVar.Inverse.InvertFor),'c')  % C
    
    if contains(lower(CtrlVar.Inverse.InvertFor),'logc')
        F.C=10.^p;
    else
        F.C=p;
    end
else
    fprintf(' CtrlVar.Inverse.InvertFor=%s \n',CtrlVar.Inverse.InvertFor)
    fprintf(' CtrlVar.Inverse.InvertFor does not have an expected value.\n')
    error('JGH:incorrect inputs')
end


[UserVar,RunInfo,F,l,drdu,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
[R,dRdp,ddRddp,RegOuts]=Regularisation(UserVar,CtrlVar,MUA,BCs,F,l,GF,Priors,Meas,BCsAdjoint,RunInfo) ;
[I,dIdp,ddIddp,MisfitOuts]=Misfit(UserVar,CtrlVar,MUA,BCs,F,l,GF,Priors,Meas,BCsAdjoint,RunInfo,drdu) ;



if nargout>1
    dJdp=dRdp+dIdp;
    Hessian=ddRddp+ddIddp;
end

if RunInfo.Forward.Converged
    ubP=F.ub;
    vbP=F.vb;
else
    ubP=[];
    vbP=[];
    I=NaN;
    MisfitOuts.I=NaN;
end

J=R+I;

if nargout>3
    JGHouts.dRdp=dRdp;
    JGHouts.dIdp=dIdp;
    JGHouts.ddIdpp=ddIddp;
    JGHouts.ddRdpp=ddRddp;
    JGHouts.RegOuts=RegOuts;
    JGHouts.MisfitOuts=MisfitOuts;
else
    JGHouts=[];
end



end

