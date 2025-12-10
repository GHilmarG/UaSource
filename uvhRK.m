 
function [Ruvh,K]=uvhRK(x,UserVar,RunInfo,CtrlVar,MUA,F0,F1)


n=MUA.Nnodes;
F1.ub=x(1:n) ;
F1.vb=x(n+1:2*n) ;
F1.h=x(2*n+1:3*n) ;

% x=[F1.ub;F1.vb;F1.h];



CtrlVar.uvhMatrixAssembly.ZeroFields=0 ; 

if nargout==1 
    CtrlVar.uvhMatrixAssembly.Ronly=1;
else
     CtrlVar.uvhMatrixAssembly.Ronly=0;
end

[UserVar,RunInfo,Ruvh,K]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1);

% Just need to normalize Ruvh

% CalcCostFunctionNRuvh(UserVar,RunInfo,CtrlVar,MUA,F1,F0,dub,dvb,dh,dl,L,luvh,cuvh,gamma,Fext0) ;

end