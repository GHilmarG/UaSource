function [Ruv,Kuv,RunInfo]=uvRK(x,CtrlVar,MUA,F,RunInfo)

nargout(3,3)

F.ub=x(1:MUA.Nnodes);
F.vb=x(MUA.Nnodes+1:2*MUA.Nnodes);


% ZeroFields=false;
% [Ruv,Kuv]=KRTFgeneralBCs(CtrlVar,MUA,F,ZeroFields) ;

CtrlVar.uvMatrixAssembly.Ronly=false;
[RunInfo,Ruv,Kuv]=uvMatrixAssembly(RunInfo,CtrlVar,MUA,F);






end