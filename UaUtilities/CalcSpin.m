function  Wxy=CalcSpin(MUA,u,v)

%
% Calculates spin around the z axis, i.e. the xy component of the spin tensor
% Wxy=CalcSpin(MUA,u,v)
%
% The spin is equal to the roation in anticlockwise direction in units rad/time. 
% If, for example, the velocity is given in the units m/yr, then Wxy has the units rad/yr.
%
% Example:
%     Wxy=CalcSpin(MUA,ub,vb);
%     figure ; PlotMeshScalarVariable(CtrlVar,MUA,Wxy*180/pi);
%     title(colorbar,'(deg./yr)') ; title('Rotation (anticlockwise) ')
% %

[~,dudy]=calcFEderivativesMUA(u,MUA);
[dvdx,~]=calcFEderivativesMUA(v,MUA);

Wxy=0.5*(dvdx-dudy);

Wxy=ProjectFintOntoNodes(MUA,Wxy);


end
