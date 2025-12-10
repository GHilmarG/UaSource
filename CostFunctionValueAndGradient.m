function [J,Gradient,Hessian,Idata,IReg,IBarrier,ub,vb,ud,vd,l,UserVar]=CostFunctionValueAndGradient(...
    UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlen,C,n,m,alpha,rho,rhow,g,GF,BCsAdjoint,Priors,Meas)
    


persistent ubP vbP PHessian

if ~isempty(ubP)
    ub=ubP;
    vb=vbP;
end


if CtrlVar.OnlyGetPersistenValues
    J=[] ; Gradient=[] ;Hessian=[] ; Idata=[] ; IReg=[] ; IBarrier=[] ; l=[];
    return
    
end


[J,Idata,IRegC,IRegAGlen,IBarrierC,IBarrierAGlen,ub,vb,ud,vd,l,dIdu,Kuv,Ruv,RunInfo]=...
    CalcMisfitFunction(UserVar,CtrlVar,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,AGlen,C,n,m,alpha,rho,rhow,g,GF,Priors,Meas);


iA=strfind(CtrlVar.AdjointGrad,'A'); iC=strfind(CtrlVar.AdjointGrad,'C'); isAgrad=~isempty(iA); isCgrad=~isempty(iC);

if isAgrad
    IReg=IRegAGlen;
    IBarrier=IBarrierAGlen;
elseif isCgrad
    IReg=IRegC;
    IBarrier=IBarrierC;
end


if nargout>1
    
    [UserVar,dJdC,dJdAGlen,ub,vb,ud,vd,xAdjoint,yAdjoint,dIdCreg,dIdAGlenreg,dIdCdata,dIdAGlendata,dIdCbarrier,dIdAGlenbarrier,lambdaAdjoint]=...
    AdjointGradientNR2d(UserVar,CtrlVar,MUA,BCs,BCsAdjoint,s,b,h,S,B,ub,vb,ud,vd,l,AGlen,C,n,m,alpha,rho,rhow,g,GF,Priors,Meas);

    
    if isAgrad
        Gradient=dJdAGlen;
        RegularisationGradient=dIdAGlenreg;
        MisfitGradient=dIdAGlendata;
    elseif isCgrad
        Gradient=dJdC;
        RegularisationGradient=dIdCreg;
        MisfitGradient=dIdCdata;
    end
    
else
    Gradient=[];
end


if nargout>2
    if isempty(PHessian)
        PHessian=MassMatrix2D1dof(MUA);
    end
    Hessian=PHessian;
else
    Hessian=[];
end

ubP=ub ; vbP=vb;




%     figure('name','FunctionGradients')
%     subplot(1,3,1) ; [pPatch,pCol]=PlotElementBasedQuantities(coordinates,connectivity,CostFunctionGradient); title('CostFunctonGradient matlab') ;
%     subplot(1,3,2) ; [pPatch,pCol]=PlotElementBasedQuantities(coordinates,connectivity,MisfitGradient); title('MisfitGradient matlab') ;
%     subplot(1,3,3) ; [pPatch,pCol]=PlotElementBasedQuantities(coordinates,connectivity,RegularisationGradient); title('RegularisationGradient matlab') ;
%
%
end

