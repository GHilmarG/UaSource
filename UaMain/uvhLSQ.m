



function [UserVar,RunInfo,F1,l1,BCs1,dt]=uvhLSQ(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1)


dt=CtrlVar.dt;
x=[F1.ub;F1.vb;F1.h];
[Aeq,beq,l]=AssembleLuvhSSTREAM(CtrlVar,MUA,BCs1,l1);

if isempty(l)
    l=beq*0 ;
end

%% find a scale for the problem

CtrlVar.uvhMatrixAssembly.ZeroFields=true ; CtrlVar.uvhMatrixAssembly.Ronly=true ; 
[UserVar,RunInfo,R0]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);   
scale=sqrt(1/(R0'*R0/2)); 

%%

fun = @(x) uvhRK(x,UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1,scale) ;


options = optimoptions('lsqnonlin','Display','iter','MaxIterations',25,'SpecifyObjectiveGradient',true,...
    'FunctionTolerance',1e-10,'Algorithm','trust-region-reflective');
%options.PlotFcn={@optimplotfirstorderoptUa,@optimplotresnormUa};
options.OptimalityTolerance = 1.000000e-15 ; options.StepTolerance = 1.000000e-20 ;
lb=[] ; ub=[] ; A=[] ; b=[] ;   nonlcon=[] ;
[x,resnorm,residual,exitflag,outputM] = lsqnonlin(fun,x,lb,ub,A,b,Aeq,beq,nonlcon,options);


n=MUA.Nnodes;
F1.ub=x(1:n) ;
F1.vb=x(n+1:2*n) ;
F1.h=x(2*n+1:3*n) ;

% [x,l,R2,r2,Qslope0,dxNorm,dlNorm,residual,g,h,output] = lsqNewton(CtrlVar,fun,x,l,Aeq,beq) ; 


end




function [Ruvh,K]=uvhRK(x,UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1,scale)


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

[UserVar,RunInfo,Ruvh,K]=uvhAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);


Ruvh=scale*Ruvh;
K=scale*K;



end




