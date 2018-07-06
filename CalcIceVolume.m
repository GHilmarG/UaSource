function [TotalIceVolume,ElementIceVolume]=CalcIceVolume(CtrlVar,MUA,h)

%  [TotalIceVolume,ElementIceVolume]=CalcIceVolume(CtrlVar,MUA,h)
%  calculates ice volume over the whole FE mesh
%  h is ice thickenss


ElementIceVolume=FEintegrate2D(CtrlVar,MUA,h) ;
TotalIceVolume=sum(ElementIceVolume);


end
