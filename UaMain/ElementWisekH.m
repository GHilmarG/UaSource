function kHe=ElementWisekH(CtrlVar,MUA,F,alpha)
%%
% kHe=ElementWisekH(CtrlVar,MUA,F)
% kHe=ElementWisekH(CtrlVar,MUA,F,alpha)
%
% (1 Nov 2026) Element-wise kH for the grounding-line mask at the integration points:
%
%   kH_e = min(CtrlVar.kH, 2/(alpha*range of h-hf over element e)),    alpha=CtrlVar.kHElementAlpha (default)
%
% where hf=rhow*max(S-B,0)/rho is the flotation thickness, evaluated from the state F. Returns a column vector of length MUA.Nele, or
% [] if alpha<=0 (the scalar CtrlVar.kH is then used everywhere).
%
% If the range of h-hf over a grounding-line element is much larger than the transition width 2/kH, the integration points practically
% never lie within the transition, the mask is effectively a step function of the thickness at the integration points, and the discrete
% uvh problem may have no root near the current iterate (the Newton iteration stagnates). kH_e widens the transition, element by element,
% to at least alpha times the range of h-hf over the element. Where the range is less than 2/(alpha*kH), kH_e=kH. The mask is only
% changed near flotation, as elsewhere it is 0 or 1 to machine precision for any kH_e.
%
% In transient runs this is evaluated once per time step from the state at the start of the time step (F0), stored in F.kHElement, and
% used by both the uvh and the uv assemblies, so that -uvh- and -uv-h- solve the same discrete problem and kH_e is fixed during the
% Newton iterations. The nodal mask (Calc_bs_From_hBS.m) keeps the scalar kH. See also ElementWisekHFromNodal.m.
%

if nargin<4 || isempty(alpha)
    alpha=0 ; if isfield(CtrlVar,"kHElementAlpha") ; alpha=CtrlVar.kHElementAlpha ; end
end

if alpha<=0
    kHe=[] ;
    return
end

d=F.h-F.rhow*max(F.S-F.B,0)./F.rho ;
kHe=ElementWisekHFromNodal(CtrlVar,reshape(d(MUA.connectivity),MUA.Nele,MUA.nod),alpha) ;

end
