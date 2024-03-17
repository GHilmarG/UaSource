

% this is just and example for how to calculate basal melt due to basal sliding
% use this as a template and change as needed


%% Get some run data 
load Ua2D_Restartfile.mat
CtrlVar=CtrlVarInRestartFile;



%%  basal melt rates due to sliding
[tbx,tby,tb,beta2] = CalcBasalTraction(CtrlVar,MUA,F.ub,F.vb,F.C,F.m,F.GF);

figure; PlotMeshScalarVariable(CtrlVar,MUA,tb) ; title(' tb (kPa)')

L=333.44 ; % Enthalpy of fusion = L = [J/gramm = kJ/kg]
ab=(tbx.*F.ub+tby.*F.vb)./(F.rho.*L);

figure ; [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,ab) ; title(' basal melt rate ')
hold on
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'r');

title(cbar,'(m/a)')
xlabel(' xps (km)') 
ylabel(' yps (km)') 


%% temperature changes due to horizontal internal deformation
%
%  dT/dt = \dot{eps} tau/(C rho)
%
%  c is here the specific heat capacity of ice
%  c=2.09 kJ/(K kg)  is the specific heat capacity of ice just below 0 C
%
% \dot{eps} tau/(c rho)   = 1/a   kJ/m^3   K  kg/kJ   m^3/kg
%                         = K/a
%
%
[etaInt,xint,yint,exx,eyy,exy,Eint,e,txx,tyy,txy]=calcStrainRatesEtaInt(CtrlVar,MUA,F.ub,F.vb,F.AGlen,F.n); 
c=2.093  ; % kJ/(K kg)


[StressPowerxy]=ProjectFintOntoNodes(MUA,exx.*txx+2*exy.*txy+eyy.*tyy);
% to include vertical directoin
ezz=-(exx+eyy) ; % in the SSA this is constant across depth 
tzz=-(txx+tyy) ; 
[StressPowerxyz]=ProjectFintOntoNodes(MUA,exx.*txx+2*exy.*txy+eyy.*tyy+ezz.*tzz);



dTdtxy=StressPowerxy./(F.rho*c);
figure ; [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,dTdtxy)  ;
hold on
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'r');
title(cbar,'(C^o/a)')
xlabel(' xps (km)')
title('dT/dt due to ice deformation in the horizontal only')
ylabel(' yps (km)') 


dTdtxyz=StressPowerxyz./(F.rho*c);
figure ; [~,cbar]=PlotMeshScalarVariable(CtrlVar,MUA,dTdtxyz)  ;
hold on
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],'r');
title(cbar,'(C^o/a)')
xlabel(' xps (km)')
title('dT/dt due to ice deformation (xyz)')
ylabel(' yps (km)') 


