function kHe=IntegrationPointkH(CtrlVar,MUA,F)
%%
% kHe=IntegrationPointkH(CtrlVar,MUA,F)
%
% (2 Nov 2026) kH=1/l for the grounding-line mask at the integration points: an array of size MUA.Nele x MUA.nip, with the nodal transition
% length F.lGL (GroundingLineTransitionLength.m) interpolated to the integration points, or, if F.lGL is empty, the scalar
% 1/CtrlVar.GroundingLineTransition.lmin. Inside a loop over integration points, use kHe(:,min(Iint,end)), which works in both cases.
%
% All routines that evaluate the grounding-line mask at integration points (the uv and uvh assemblies, and the gradient and Hessian
% calculations of the inverse problem) use this, so that they are consistent. The nodal mask (Calc_bs_From_hBS.m) and the water-depth mask
% use the scalar 1/lmin.
%

kHe=1/CtrlVar.GroundingLineTransition.lmin ;
lGL=[] ; if (isobject(F) && isprop(F,"lGL")) || (isstruct(F) && isfield(F,"lGL")) ; lGL=F.lGL ; end
if numel(lGL)~=MUA.Nnodes
    return
end
lnod=reshape(lGL(MUA.connectivity),MUA.Nele,MUA.nod) ;
kHe=zeros(MUA.Nele,MUA.nip) ;
for Iint=1:MUA.nip
    fun=shape_fun(Iint,2,MUA.nod,MUA.points) ;
    kHe(:,Iint)=1./(lnod*fun) ;
end

end
