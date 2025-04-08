function gphi=DegradationFunction(CtrlVar,phi)

%  phi=0 : undamaged 
%  phi=1 : fully damaged 
%
%  gphi=0 : fully damaged
%  gphi=1 : undamaged
%
narginchk(2,2)



k=CtrlVar.PhaseFieldFracture.k ;

gphi=(1-k)* (1-phi).^2 + k ;


end

