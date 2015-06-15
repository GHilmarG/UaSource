function [J,CostFunctionGradient,Idata,IReg,IBarrier,ub,vb,ud,vd,lx,ly] = ...
    CostFunctionValueAndGradient(CtrlVar,MUA,JoptVector,...
    s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,...
    sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,...
    AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
    Luv,Luvrhs,lambdauv,LAdjoint,LAdjointrhs,lambdaAdjoint,...
    GF)


[J,Idata,IRegC,IRegAGlen,IBarrierC,IBarrierAGlen,ub,vb,ud,vd,dIdu,K]=...
    CalcMisfitFunction(CtrlVar,MUA,s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,...
    sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,...
    AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
    Luv,Luvrhs,lambdauv,GF);



iA=strfind(CtrlVar.AdjointGrad,'A'); iC=strfind(CtrlVar.AdjointGrad,'C'); isAgrad=~isempty(iA); isCgrad=~isempty(iC);

if isAgrad
    IReg=IRegAGlen;
    IBarrier=IBarrierAGlen;
elseif isCgrad
    IReg=IRegC;
    IBarrier=IBarrierC;
end


if nargout>1
    
    [dJdC,dJdAGlen,ub,vb,ud,vd,lx,ly,dIdCreg,dIdAGlenreg,dIdCdata,dIdAGlendata,dIdCbarrier,dIdAGlenbarrier]=...
        AdjointGradientNR2d(CtrlVar,MUA,s,b,h,S,B,ub,vb,ud,vd,AGlen,n,C,m,rho,rhow,alpha,g,...
        sMeas,uMeas,vMeas,wMeas,bMeas,BMeas,...
        AGlen_prior,CAGlen,C_prior,CC,b_prior,Cd,...
        Luv,Luvrhs,lambdauv,LAdjoint,LAdjointrhs,lambdaAdjoint,GF,K);
    
    
    if isAgrad
        CostFunctionGradient=dJdAGlen;
        RegularisationGradient=dIdAGlenreg;
        MisfitGradient=dIdAGlendata;
    elseif isCgrad
        CostFunctionGradient=dJdC;
        RegularisationGradient=dIdCreg;
        MisfitGradient=dIdCdata;
    end
    
else
    CostFunctionGradient=[];
end

%     figure('name','FunctionGradients')
%     subplot(1,3,1) ; [pPatch,pCol]=PlotElementBasedQuantities(coordinates,connectivity,CostFunctionGradient); title('CostFunctonGradient matlab') ;
%     subplot(1,3,2) ; [pPatch,pCol]=PlotElementBasedQuantities(coordinates,connectivity,MisfitGradient); title('MisfitGradient matlab') ;
%     subplot(1,3,3) ; [pPatch,pCol]=PlotElementBasedQuantities(coordinates,connectivity,RegularisationGradient); title('RegularisationGradient matlab') ;
%
%
end

