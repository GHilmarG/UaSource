function gphi=DegradationFunction(CtrlVar,phi)

% $\phi$ is the phase-field variable.
%
% $\phi=0$ for undamaged material 
% $\phi=1$ for (fully) damaged material 
%
%
%
% $$ g_{\phi}=(1-k) (1-\phi)^2 + k $$
%
%
%  gphi=0 : fully damaged
%  gphi=1 : undamaged
%

narginchk(2,2)



k=CtrlVar.PhaseFieldFracture.k ;

gphi=(1-k)* (1-phi).^2 + k ;


end

