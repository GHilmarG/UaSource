



function [UserVar,RunInfo,F1,l1,BCs1,dt]=uvhLSQ(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1)



x=[F1.ub;F1.vb;F1.h];
[Aeq,beq,luvh0]=AssembleLuvhSSTREAM(CtrlVar,MUA,BCs1,l1);

if isempty(luvh0)
    l1=cuvh*0 ;
end

fun = @(x) uvhRK(x,UserVar,RunInfo,CtrlVar,MUA,F0,F1) ;


[x,l,R2,r2,Qslope0,dxNorm,dlNorm,residual,g,h,output] = lsqNewton(CtrlVar,fun,x,l1,Aeq,beq) ; 


end




function [Ruvh,K]=uvhRK(x,UserVar,RunInfo,CtrlVar,MUA,F0,F1)


n=MUA.Nnodes;
F1.ub=x(1:n) ;
F1.vb=x(n+1:2*n) ;
F1.h=x(2*n+1:3*n) ;


CtrlVar.uvhMatrixAssembly.ZeroFields=0 ; 

if nargout==1 
    CtrlVar.uvhMatrixAssembly.Ronly=1;
else
    CtrlVar.uvhMatrixAssembly.Ronly=0;
end

[UserVar,RunInfo,Ruvh,K]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1);



end




