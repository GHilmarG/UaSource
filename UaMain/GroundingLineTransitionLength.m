function lGL=GroundingLineTransitionLength(CtrlVar,MUA,F,dt)
%%
% lGL=GroundingLineTransitionLength(CtrlVar,MUA,F)
% lGL=GroundingLineTransitionLength(CtrlVar,MUA,F,dt)
%
% (2 Nov 2026) Nodal transition length of the grounding-line mask, G=(1+tanh((h-hf)/l))/2 at the integration points (l interpolated from
% the nodes), in units of thickness (m):
%
%   l = min(lmax, max(lmin, c*dt*|dh/dt|))
%
% with lmin, lmax and c from CtrlVar.GroundingLineTransition, dh/dt=F.dhdt the rate of thickness change of the previous time step, and dt the
% current time step (default CtrlVar.dt).
%
% Why: the uvh Newton iteration converges quadratically only if the change of h-hf over the time step at a point stays within the range
% over which the mask can be linearised, ie is small compared with the transition length. With a fixed l=lmin (1 m), the iteration
% stagnated when the grounding line moved by more than a small fraction of an element per time step: the mask was then effectively a step
% function of the thickness at the integration points, and the discrete problem could have no root near the iterate. Making l a fixed
% fraction c of the expected thickness change over the time step removes this limitation on the time step. The smoothing is proportional
% to dt (like the error of the implicit time stepping) and vanishes as dt->0. lmax limits the smoothing where |dh/dt| is large far from
% flotation (eg near calving fronts), where a large l would give a spurious non-zero mask on floating ice.
%
% l is evaluated once per time step from the state at its start (uvhRootFinding.m, uvhSemiImplicit.m; recalculated if the time step is
% reduced, uvh2NotConvergent.m), stored in F.lGL, and kept fixed during the Newton iterations, so the Jacobian has no additional terms.
% It is used by the uvh and uv assemblies (vectorised and compiled) and by the inverse problem, through IntegrationPointkH.m.
%
% Returns [] (ie the scalar lmin everywhere) if c=0, if there is no positive time step or no rate of thickness change, or if l=lmin at all
% nodes.
%

GLT=CtrlVar.GroundingLineTransition ;
if nargin<4 || isempty(dt) ; dt=CtrlVar.dt ; end
lGL=[] ;
rate=[] ; if (isobject(F) && isprop(F,"dhdt")) || (isstruct(F) && isfield(F,"dhdt")) ; rate=F.dhdt ; end
if GLT.c<=0 || isempty(dt) || ~(dt>0) || numel(rate)~=MUA.Nnodes
    return
end
rate(~isfinite(rate))=0 ;
l=min(GLT.lmax,max(GLT.lmin,GLT.c*dt*abs(rate(:)))) ;
if all(l==GLT.lmin)
    return
end
lGL=l ;
if CtrlVar.InfoLevel>=10
    fprintf(" Grounding-line transition length: median %.3g m, max %.3g m, %i nodes at lmax=%g m \n",median(lGL),max(lGL),nnz(lGL>=GLT.lmax),GLT.lmax) ;
end

end
