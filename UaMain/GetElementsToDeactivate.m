



function [UserVar,ElementsToBeDeactivated]=GetElementsToDeactivate(UserVar,RunInfo,CtrlVar,MUA,F,BCs,ElementsToBeDeactivated)


narginchk(7,7)
nargoutchk(2,2)

InputFile="DefineElementsToDeactivate.m";
TestIfInputFileInWorkingDirectory(InputFile) ;

NarginInputFile=nargin(InputFile);
NargoutInputFile=nargout(InputFile);


if NarginInputFile==18

% This is the old input format, which is still accepted

    [UserVar,ElementsToBeDeactivated]=...
        DefineElementsToDeactivate(UserVar,RunInfo,CtrlVar,MUA,MUA.xEle,MUA.yEle,ElementsToBeDeactivated,F.s,F.b,F.S,F.B,F.rho,F.rhow,F.ub,F.vb,F.ud,F.vd,F.GF);

elseif NarginInputFile==7

    [UserVar,ElementsToBeDeactivated]=DefineElementsToDeactivate(UserVar,RunInfo,CtrlVar,MUA,F,BCs,ElementsToBeDeactivated);


else

    fprintf("The function declaration of DefineElementsToDeactive.m must be on the form: \n")
    fprintf("    [UserVar,ElementsToBeDeactivated]=DefineElementsToDeactivate(UserVar,RunInfo,CtrlVar,MUA,F,bcS,ElementsToBeDeactivated) \n")
    error('Ua:GetElementsToDeactive','DefineElementsToDeactive must have 7 input  arguments')

end





end