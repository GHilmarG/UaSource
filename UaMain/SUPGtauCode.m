function Code=SUPGtauCode(tauOption)

%% Integer code of the SUPG tau option (CtrlVar.uvh.SUPG.tau), as used by SUPGtau.m
%
%   Code=SUPGtauCode(CtrlVar.uvh.SUPG.tau)
%
%   Code   Option
%    1     "tau1"   coth(ECN)-1/ECN times the spatial definition (textbook recommendation for linear advection)
%    2     "tau2"   inversely weighted average of the spatial and the temporal definition
%    3     "taus"   spatial definition, l/(2v), independent of the time step
%    4     "taut"   temporal definition, dt/2, independent of the speed
%    5     "tau0"   zero
%
% SUPGtau.m accepts the option either as the name or as this code. The compiled assembly (UaSource/UaCompiled) determines the code
% once per assembly (MakeKernelParamsFlex.m), instead of comparing strings at every integration point.
%
% Written so that it can be used by MATLAB Coder (character comparisons only).
%
% See also: SUPGtau
%
% (26 Oct 2026)
%%

L=char(tauOption);

if strcmp(L,'tau1')
    Code=1;
elseif strcmp(L,'tau2')
    Code=2;
elseif strcmp(L,'taus')
    Code=3;
elseif strcmp(L,'taut')
    Code=4;
elseif strcmp(L,'tau0')
    Code=5;
else
    error('SUPGtauCode:UnknownOption','Unknown SUPG tau option "%s" (CtrlVar.uvh.SUPG.tau). Options: "tau1", "tau2", "taus", "taut", "tau0".',L)
end

end
