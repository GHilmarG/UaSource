function [tbx,tby,tb] = CalcBasalTraction(CtrlVar,UserVar,MUA,F)

narginchk(4,4)

 % [tbx,tby,tb,beta2] = CalcBasalTraction(CtrlVar,MUA,ub,vb,C,m,GF)  ; % old
 % version

%%
%
%    [tbx,tby,tb,beta2] = CalcBasalTraction(CtrlVar,MUA,ub,vb,C,m,GF)
% 
% Calculates basal traction from basal velocity using the sliding law.
% 
% Returns nodal values 
%
% Note: This can only be used to calculate basal traction when using the SSTREAM
% and the Hybrid flow approximation. This will not return correct results for
% the SSHEET approximation!
%
% Note: There is a slight inconsistency with respect to how this is done
% internally in Ua in the sense that the floating mask is here evalutated at
% nodes, whereas internally this is done at integration points. 
%
%
%%


% 
% if CtrlVar.CisElementBased
%     % project onto nodes
%     [M,ElePerNode] = Ele2Nodes(MUA.connectivity,MUA.Nnodes); 
% 
%     F.C=M*F.C;
%     F.m=M*F.m;
% 
% end
% 
% beta2=(F.C+CtrlVar.Czero).^(-1./F.m).*(sqrt(F.ub.*F.ub+F.vb.*F.vb+CtrlVar.SpeedZero^2)).^(1./F.m-1) ;
%
% tbx=F.GF.node.*beta2.*F.ub;
% tby=F.GF.node.*beta2.*F.vb;
% tb=sqrt(tbx.^2+tby.^2);


hf=F.rhow*(F.S-F.B)./F.rho ;
He = HeavisideApprox(CtrlVar.kH,F.h-hf,CtrlVar.Hh0);  % 1
delta = DiracDelta(CtrlVar.kH,F.h-hf,CtrlVar.Hh0) ;

[tbx,tby] = ...
    BasalDrag(CtrlVar,MUA,He,delta,F.h,F.B,F.S-F.B,F.rho,F.rhow,F.ub,F.vb,F.C,F.m,F.uo,F.vo,F.Co,F.mo,F.ua,F.va,F.Ca,F.ma,F.q,F.g,F.muk);

tb=sqrt(tbx.^2+tby.^2);


end

