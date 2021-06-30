function [VAF,IceVolume,GroundedArea,hAF,hfPos]=CalcVAF(CtrlVar,MUA,h,B,S,rho,rhow,GF)

%%
%
%   [VAF,IceVolume,GroundedArea,hAF,hfPos]=CalcVAF(CtrlVar,MUA,h,B,S,rho,rhow,GF)
%
% Calculates volume above flotation, and optionally ice volume and grounded area
%
% GF is only needed to calculate grounded area.
%
%
% To calculate a rough estimate of resulting change in mean sea level, divide the change
% in VAF with the area of the ocean (3.625e14 m^2).
%
%
%   VAF       :  Volume above flotation
%   IceVolume :  Total ice volume withing the domain, i.e. including areas that are afloat.
%   hAF       :  (postive) ice thickness above floation
%   hfPOs     :  (positive) flotation thickness (also somtimes referred to as floation profile). Where h>fhPos, the ice is grounded.
%
% Example:
%
%   load("PIG-TWG-RestartFile.mat") 
%   [VAF,IceVolume,GroundedArea,hAF,hfPos]=CalcVAF([],MUA,F.h,F.B,F.S,F.rho,F.rhow,F.GF);
%   CtrlVar=CtrlVarInRestartFile;
%   FindOrCreateFigure("VAF") ; 
%   [~,cbar]=PlotMeshScalarVariable(CtrlVarInRestartFile,MUA,hAF) ; 
%   axis tight
%   hold on ; PlotLatLonGrid(CtrlVar.PlotXYscale) ;
%   hold on ; PlotGroundingLines(CtrlVar,MUA,F.GF,[],[],[],'r');
%   xlabel("xps (km)") ; ylabel("yps (km)") ; title(cbar,"(m)") ; title("ice thickness above flotation")
%   fprintf("VAF=%f (Gt/yr)\n",VAF.Total/1e9)   ; 
%   fprintf("GroundedArea=%-7.2f (times the area of iceland)\n",GroundedArea.Total/1e6/103e3) ; 
%%

narginchk(7,8)
nargoutchk(1,5)


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
hfPos=(S>B).*rhow.*(S-B)./rho ;            % (positive) flotation thickness
hAF= (h>hfPos).*(h-hfPos) ;                % (positive) ice thickness above floatation



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