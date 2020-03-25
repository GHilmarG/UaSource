function [VAF,IceVolume,GroundedArea]=CalcVAF(CtrlVar,MUA,h,B,S,rho,rhow,GF)

%%
% Calculates volume above flotation, and optionally ice volume and grounded area
%
% GF is only needed to calculate volume and grounded area. 
%
%
% To calculate a rough estimate of resulting change in mean sea level, divide the change
% in VAF with the area of the ocean (3.625e14 m^2).
%%

narginchk(7,8)
nargoutchk(1,3)


% One option: 
%
% hf=rhow*(S-B)./rho ;
% hAF=h*0 ;                                   % hAF : ice-thickness above flotation. 
% 
% isBgtS= B > S ;  hAF(isBgtS)=h(isBgtS) ;    % grounded above sea level, full contribution to hAF
% I=~isBgtS & ~ishlthf ; hAF(I)=h(I)-hf(I) ;  % grounded below sea level some contribution to hAF
% ishlthf = (h < hf)  ; hAF(ishlthf)=0   ;    % not grounded, no contribution to hAF
% 
% 
% or simply: 
hfPos=(S>B).*rhow.*(S-B)./rho ; % (positive) flotation thickness
hAF= (h-hfPos) ;                % ice thickness above floatation


VAF.node=hAF.*rho./rhow ;               % thickness above flotation in water eq. 
VAF.ele=FEintegrate2D([],MUA,VAF.node); % VAF for each element (m^3) 
VAF.Total=sum(VAF.ele);                 % total volume above flotation over the whole model domain 


if nargout>1
    IceVolume.Ele=FEintegrate2D([],MUA,h);
    IceVolume.Total=sum(IceVolume.Ele);
    
    GroundedArea.Ele=FEintegrate2D([],MUA,GF.node);
    GroundedArea.Total=sum(GroundedArea.Ele);
end




end