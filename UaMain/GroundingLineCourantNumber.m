function [CoMedian,CoP90,CoMax,nGL]=GroundingLineCourantNumber(CtrlVar,MUA,F0,F1)
%%
% [CoMedian,CoP90,CoMax,nGL]=GroundingLineCourantNumber(CtrlVar,MUA,F0,F1)
%
% (2 Nov 2026) Grounding-line Courant number of a time step: for each element crossing flotation at the start (F0) or at the end (F1) of
% the time step,
%
%   Co_e = max over the nodes of the element of |(h-hf) at the end - (h-hf) at the start| / range of (h-hf) over the element at the start,
%
% where hf=rhow*max(S-B,0)/rho is the flotation thickness. Since the range of h-hf over an element is about |grad(h-hf)| times the element
% size, and the change of h-hf at a point is about |grad(h-hf)| times the grounding-line displacement, Co_e is approximately the
% grounding-line displacement over the time step in units of the element size (the grounding-line Courant number).
%
% It is a diagnostic of how fast the grounding line moves relative to the mesh. With a fixed transition length of the grounding-line mask
% (proportional to the range of h-hf over an element), the uvh Newton iteration converged reliably only up to Courant numbers of a few
% percent; with the nodal transition length proportional to the thickness change over the time step (GroundingLineTransitionLength.m)
% this limitation is removed.
%
% Returns the median, the 90th percentile and the maximum of Co_e over these elements, and their number nGL (NaN if there are none). F0 and
% F1 must be on the same mesh. For elements with more than three nodes the nodal range is used as an approximation.
%

CoMedian=NaN ; CoP90=NaN ; CoMax=NaN ; nGL=0 ;
if isempty(F0.h) || isempty(F1.h) || numel(F0.h)~=MUA.Nnodes || numel(F1.h)~=MUA.Nnodes
    return
end

d0=F0.h-F0.rhow*max(F0.S-F0.B,0)./F0.rho ;
d1=F1.h-F1.rhow*max(F1.S-F1.B,0)./F1.rho ;
dn0=reshape(d0(MUA.connectivity),MUA.Nele,MUA.nod) ;
dn1=reshape(d1(MUA.connectivity),MUA.Nele,MUA.nod) ;
GL=(any(dn0>0,2) & any(dn0<0,2)) | (any(dn1>0,2) & any(dn1<0,2)) ;
nGL=nnz(GL) ;
if nGL==0
    return
end

span=max(dn0(GL,:),[],2)-min(dn0(GL,:),[],2) ;
Co=max(abs(dn1(GL,:)-dn0(GL,:)),[],2)./max(span,eps) ;
CoMedian=median(Co) ; CoP90=prctile(Co,90) ; CoMax=max(Co) ;

if CtrlVar.InfoLevel>=1
    fprintf(" Grounding-line Courant number (over %i elements): median %.3g, 90%% %.3g, max %.3g \n",nGL,CoMedian,CoP90,CoMax) ;
end

end
