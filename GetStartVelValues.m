function [UserVar,F]=GetStartVelValues(UserVar,CtrlVar,MUA,BCs,F)

narginchk(5,5)
nargoutchk(2,2)


[UserVar,F.ub,F.vb,F.ud,F.vd]=DefineStartVelValues(UserVar,CtrlVar,MUA,BCs,F.ub,F.vb,F.ud,F.vd,CtrlVar.time,F.s,F.b,F.h,F.S,F.B,F.rho,F.rhow,F.GF,F.AGlen,F.n,F.C,F.m);


end

