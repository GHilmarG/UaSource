function [VAF,IceVolume,GroundedArea]=CalcVAF(CtrlVar,MUA,h,b,S,rho,rhow,GF)

%%
% Calculates volume above flotation, and optionally ice volume and grounded area
%
% GF is only needed to calculate volume and grounded area. 
%
%
%%

narginchk(7,8)
nargoutchk(1,3)

VAF.node=(h+rhow.*(b-S)./rho).*heaviside(S-b);
VAF.ele=FEintegrate2D([],MUA,VAF.node);
VAF.Total=sum(VAF.ele);


if nargout>1
    IceVolume.Ele=FEintegrate2D([],MUA,h);
    IceVolume.Total=sum(IceVolume.Ele);
    
    GroundedArea.Ele=FEintegrate2D([],MUA,GF.node);
    GroundedArea.Total=sum(GroundedArea.Ele);
end



end