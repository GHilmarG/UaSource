function [UserVar,EleSizeDesired,ElementsToBeRefined]=GetDesiredEleSize(UserVar,CtrlVar,MUA,F,GF,x,y,EleSizeDesired,ElementsToBeRefined,NodalErrorIndicators)

narginchk(10,10)
nargoutchk(3,3)


[UserVar,EleSizeDesired,ElementsToBeRefined]=...
    DefineDesiredEleSize(UserVar,CtrlVar,MUA,x,y,EleSizeDesired,ElementsToBeRefined,F.s,F.b,F.S,F.B,F.rho,F.rhow,F.ub,F.vb,F.ud,F.vd,GF,NodalErrorIndicators);



end

