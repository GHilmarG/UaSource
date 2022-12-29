


function [VAF,IceVolume,GroundedArea,hAF,hfPos]=CalcVAF(CtrlVar,MUA,h,B,S,rho,rhoOcean,GF,options)

%%
%
%   [VAF,IceVolume,GroundedArea,hAF,hfPos]=CalcVAF(CtrlVar,MUA,h,B,S,rho,rhow,GF)
%
%
%   rhoOcean is here the density of the "water" in which the ice is floating.
%   This would most likely be the density of the ocean
%
% Calculates volume above flotation, and optionally ice volume and grounded area
%
%  VAF has the units distance^3, i.e. it is a volume, not weight.  It is the water equivalent volume.
%
% If all distance units are in meters, and we divide VAF as calculated by 10e9, then the units of VAF are km^3
%
%
% GF is only needed to calculate grounded area.
%
%
% To calculate a rough estimate of resulting change in mean sea level, divide the change in VAF with the area of the ocean
% (3.625e14 m^2). Since 1Gt is = 1e9 m^3 water equivalent, the conversion between sea-level change and ice loss is about
%
%     0.001/362.5  (m/Gt)
%
% so about 1 mm sea level change for every 362.5 Gt water added.  This is the sea level potential per Gt water. 
% 
% This calculation does not account for other effecs such as
% ocean salinity changes, but these are only expected to change the value by a few %.
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
%   [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,hAF) ; 
%   axis tight
%   hold on ; PlotLatLonGrid(CtrlVar.PlotXYscale) ;
%   hold on ; PlotGroundingLines(CtrlVar,MUA,F.GF,[],[],[],'r');
%   xlabel("xps (km)",interpreter="latex") ; ylabel("yps (km)",interpreter="latex") ; 
%   title(cbar,"(m)") ; title("ice thickness above flotation")
%   fprintf("VAF=%f (Gt/yr)\n",VAF.Total/1e9)   ; 
%   fprintf("GroundedArea=%-7.2f (times the area of iceland)\n",GroundedArea.Total/1e6/103e3) ; 
%   colormap(othercolor('Blues7',1024));
%
%%


arguments
    CtrlVar     struct
    MUA         struct
    h           (:,1)  double
    B           (:,1)  double
    S           (:,1)  double
    rho         (:,1)  double
    rhoOcean    (:,1)  double
    GF          struct
    options.boundary double=nan
    options.plot logical = false 
end


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
hfPos=(S>B).*rhoOcean.*(S-B)./rho ;            % (positive) flotation thickness
hAF= (h>hfPos).*(h-hfPos) ;                    % (positive) ice thickness above floatation

if ~isnan(options.boundary)  % OK boundary was given as input, so only calculate VAF inside of that boundary
    xy=[MUA.coordinates(:,1) MUA.coordinates(:,2)] ;
    isInside=inpoly2(xy,options.boundary);
    hAF(~isInside)=0;                       % simply set all nodal values outside of that boundary to zero. 
end

VAF.node=hAF.*rho./rhoOcean ;                % thickness above flotation in (ocean) water equivalent.


VAF.ele=FEintegrate2D(CtrlVar,MUA,VAF.node); % VAF for each element (m^3)
VAF.Total=sum(VAF.ele);                      % total volume above flotation over the whole model domain




if nargout>1
    IceVolume.Ele=FEintegrate2D(CtrlVar,MUA,h);
    IceVolume.Total=sum(IceVolume.Ele) ;  % ice voluem volume
                                          % to get water volume multiply with rhoice/rhowater
                                                  
    
    GroundedArea.Ele=FEintegrate2D(CtrlVar,MUA,GF.node);
    GroundedArea.Total=sum(GroundedArea.Ele);
end


if options.plot 

  FindOrCreateFigure("VAF") ; 
  [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,hAF) ; 
  axis tight
  hold on ; PlotLatLonGrid(CtrlVar.PlotXYscale) ;
  hold on ; PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'r');
  xlabel("xps (km)",interpreter="latex") ; ylabel("yps (km)",interpreter="latex") ; 
  title(cbar,"(m)") ; title("ice thickness above flotation")
  fprintf("VAF=%f (Gt/yr)\n",VAF.Total/1e9)   ; 
  fprintf("GroundedArea=%-7.2f (times the area of iceland)\n",GroundedArea.Total/1e6/103e3) ; 
  colormap(othercolor('Blues7',1024));

end


end