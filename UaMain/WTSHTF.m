





function [UserVar,RunInfo,F,F0,l]= WTSHTF(UserVar,RunInfo,CtrlVar,MUA,BCs,F,F0,Fm1,l)

%%
%
% If the implicit uvh solution did not converge, take a semi-implicit step and
% continue.
%
%  n adapt time step the dt will be reduced by the fraction
%
%     CtrlVar.ATStimeStepFactorDownNOuvhConvergence
%
% So this should hopefully cause a graceful reduction in time step without
% significant reduction in accuracy.
%

narginchk(9,9)
nargoutchk(5,5)

if ~isfield(CtrlVar,"Try_uv_SolveIf_uvh_SolveNotConvergent")
    CtrlVar.Try_uv_SolveIf_uvh_SolveNotConvergent=false; 
end


if CtrlVar.Try_uv_SolveIf_uvh_SolveNotConvergent
    fprintf("WTSHTF: Now trying uv-h solve. \n")
    [UserVar,RunInfo,F,F0,l]= uvhSemiImplicit(UserVar,RunInfo,CtrlVar,MUA,F0,F0,l,BCs)  ;
end
   

end


