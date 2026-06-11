







function [Ruvh,K]=uvhRK(x,UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1,scale)

persistent nCalls

if isempty(nCalls)
    nCalls.rhs=0;
    nCalls.Assembly=0;
end


n=MUA.Nnodes;
F1.ub=x(1:n) ;
F1.vb=x(n+1:2*n) ;
F1.h=x(2*n+1:3*n) ;


CtrlVar.uvhMatrixAssembly.ZeroFields=0 ;

if nargout==1
    CtrlVar.uvhMatrixAssembly.Ronly=1;
    nCalls.rhs=nCalls.rhs+1;
else
    CtrlVar.uvhMatrixAssembly.Ronly=0;
    nCalls.rhs=nCalls.rhs+1;
    nCalls.Assembly=nCalls.Assembly+1;

end

%fprintf("calls: rhs %i \t assembly %i \n %",nCalls.rhs,nCalls.Assembly)

[UserVar,RunInfo,Ruvh,K]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);


Ruvh=scale*Ruvh;
K=scale*K;



end

