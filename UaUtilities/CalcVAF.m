function [VAF,IceVolume,GroundedArea]=CalcVAF(CtrlVar,MUA,h,B,S,rho,rhow,GF)

%%
% Calculates volume above flotation, and optionally ice volume and grounded area
%
% GF is only needed to calculate volume and grounded area. 
%
%
%%

narginchk(7,8)
nargoutchk(1,3)

hf=(S-B).*rhow./rho ;       % if h> hf, then ice is grounded  (but note that hf is negative for B>S)
hf(hf<0)=0;                 % this is the positive flotation thickness 

VAF.node=(h-hf).*GF.node ;  % thickness above flotation over grounded nodes
% VAF.node=(h-hf).*heaviside(h-hf) 

VAF.ele=FEintegrate2D([],MUA,VAF.node);
VAF.Total=sum(VAF.ele);


if nargout>1
    IceVolume.Ele=FEintegrate2D([],MUA,h);
    IceVolume.Total=sum(IceVolume.Ele);
    
    GroundedArea.Ele=FEintegrate2D([],MUA,GF.node);
    GroundedArea.Total=sum(GroundedArea.Ele);
end



end