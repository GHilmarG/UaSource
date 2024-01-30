function [Ruv,Kuv,Tint,Fext]=KRTFgeneralBCs(CtrlVar,MUA,F)

%%
%
%   [Ruv,Kuv,Tint,Fext]=KRTFgeneralBCs(CtrlVar,MUA,F,ZeroFields)
%
%  uv SSA/SSTREAM assembly
%
%
%
%%

narginchk(3,3)



% Make sure that the info about zero fields is always in the CtrlVar on input, and not implied by the number of input
% arguments.

% if nargin<5
%     CtrlVar.uvAssembly.ZeroFields=false;
% else
%     CtrlVar.uvAssembly.ZeroFields=ZeroFields;
% end



if nargout==1 && ~CtrlVar.uvMatrixAssembly.Ronly
    error("KRTFgeneralBCs: More than one output, but assembly only required for R. ")
end


[Ruv,Kuv,Tint,Fext]=uvMatrixAssembly(CtrlVar,MUA,F);




end





