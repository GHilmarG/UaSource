function [UserVar,F,l]=GetStartVelValues(UserVar,CtrlVar,MUA,BCs,F,l)

narginchk(6,6)
nargoutchk(3,3)

InputFile="DefineStartVelValues.m" ; 

% TestIfInputFileInWorkingDirectory(InputFile) ;



N=nargout(InputFile);
NargInputFile=nargin(InputFile);

if NargInputFile>6

    [UserVar,F.ub,F.vb,F.ud,F.vd]=DefineStartVelValues(UserVar,CtrlVar,MUA,BCs,F.ub,F.vb,F.ud,F.vd,CtrlVar.time,F.s,F.b,F.h,F.S,F.B,F.rho,F.rhow,F.GF,F.AGlen,F.n,F.C,F.m);

else

    [UserVar,F.ub,F.vb,F.ud,F.vd,l]=DefineStartVelValues(UserVar,CtrlVar,MUA,BCs,F,l) ;


end

