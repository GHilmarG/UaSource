function kHe=ElementWisekHFromNodal(CtrlVar,dnod,alpha)
%%
% kHe=ElementWisekHFromNodal(CtrlVar,dnod,alpha)
%
% (1 Nov 2026) Element-wise kH from the nodal values of h-hf of each element, dnod (Nele x nod):
%
%   kH_e = min(CtrlVar.kH, 2/(alpha*(max(dnod,[],2)-min(dnod,[],2))))
%
% This is the single definition of the element-wise kH, used by ElementWisekH.m and by the assembly routines that work with nodal
% values per element. See ElementWisekH.m.
%

span=max(dnod,[],2)-min(dnod,[],2) ;
kHe=min(CtrlVar.kH,2./(alpha*max(span,1e-12))) ;

end
