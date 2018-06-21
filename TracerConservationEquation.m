function [UserVar,c1,lambda]=TracerConservationEquation(UserVar,CtrlVar,MUA,dt,c0,u0,v0,a0,u1,v1,a1,kappa,BCsTracer)





MLC=BCs2MLC(MUA,BCsTracer);
L=MLC.hL ; Lrhs=MLC.hRhs ; lambda=Lrhs*0;

[UserVar,kv,rh]=TracerConservationEquationAssembly(UserVar,CtrlVar,MUA,dt,c0,u0,v0,a0,u1,v1,a1,kappa);

[c1,lambda]=solveKApe(kv,L,rh,Lrhs,c0,lambda,CtrlVar);
c1=full(c1);


end




