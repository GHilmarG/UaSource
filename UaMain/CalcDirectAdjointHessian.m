

function  [Hsparse,Hfull,g0,J0] = CalcDirectAdjointHessian(func,p,CtrlVar,iRange)


% [J,dJdp,Hessian,JGHouts,F,RunInfo]=JGH(p,plb,pub,UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo)


[J0,g0,~,JGHouts]=func(p0);

ddRddp=JGHouts.ddRddp;
uAdjoint=JGHouts.MisfitOuts.uAdjoint;
vAdjoint=JGHouts.MisfitOuts.vAdjoint;


%{
  JGHouts.dRdp=dRdp;
    JGHouts.dIdp=dIdp;
    JGHouts.ddIdpp=ddIddp;
    JGHouts.ddRdpp=ddRddp;
    JGHouts.RegOuts=RegOuts;
    JGHouts.MisfitOuts=MisfitOuts;
%}



[dudA,dvdA]=duvdAFunc(CtrlVar,MUA,F,BCs) ;
[dudC,dvdC]=duvdCFunc(CtrlVar,MUA,F,BCs) ;




end