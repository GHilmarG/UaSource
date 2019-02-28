function [UserVar,dhdt]=Calculate_dhdt(UserVar,CtrlVar,MUA,F,BCs)
%%
% Calculates dh/dt from flux divergence as
%
%   dh/dt = a -  ( dqx/dx + dqy/dy)
%
% where qx=uh and qy =vh.
%
%   [UserVar,dhdt]=Calculate_dhdt(UserVar,CtrlVar,MUA,F,BCs)
%
% uses u=F.ub, and hence only correct for plug flow, e.g. SSA
%
% Projects the values directly onto nodes.
%
% Optionally uses SUPG, but very doubtfull that this SUPG treatment is required.
%
%

narginchk(5,5)

if nargin<5
    BCs=[];
end

if ~isfield(CtrlVar,'dhdtCalculationUsesSUPG')
    CtrlVar.dhdtCalculationUsesSUPG=false;
end

if CtrlVar.dhdtCalculationUsesSUPG
    [UserVar,dhdt]=dhdtExplicitSUPG(UserVar,CtrlVar,MUA,F,BCs);
else
    [UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs) ;
    
end

end
