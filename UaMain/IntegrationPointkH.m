function kHe=IntegrationPointkH(CtrlVar,MUA,F)
%%
% kHe=IntegrationPointkH(CtrlVar,MUA,F)
%
% (1 Nov 2026) kH for the grounding-line mask at the integration points: a column vector of length MUA.Nele (element-wise kH, see
% ElementWisekH.m) if CtrlVar.kHElementAlpha>0, otherwise the scalar CtrlVar.kH. If F.kHElement is set (in transient runs at the start of
% each time step from the state at its start, and in inversions once at the start of the inversion), it is used; otherwise kH_e is
% calculated from the state F.
%
% All routines that evaluate the grounding-line mask at integration points (the uv and uvh assemblies, and the gradient and Hessian
% calculations of the inverse problem) use this, so that they are consistent. The nodal mask (Calc_bs_From_hBS.m) keeps the scalar kH.
%

alpha=0 ; if isfield(CtrlVar,"kHElementAlpha") ; alpha=CtrlVar.kHElementAlpha ; end
if alpha<=0
    kHe=CtrlVar.kH ;
    return
end

kHeField=[] ;
if (isobject(F) && isprop(F,"kHElement")) || (isstruct(F) && isfield(F,"kHElement")) ; kHeField=F.kHElement ; end
if numel(kHeField)==MUA.Nele
    kHe=kHeField(:) ;
else
    kHe=ElementWisekH(CtrlVar,MUA,F,alpha) ;
end

end
