





function Fnew=PFFmapF(UserVar,RunInfo,CtrlVar,MUAold,MUAnew,F)

Fnew=F ;


Fnew.x=MUAnew.coordinates(:,1);
Fnew.y=MUAnew.coordinates(:,2);


Fnew.AGlen=F.AGlen(1)+zeros(MUAnew.Nnodes,1);
Fnew.AGlen0=F.AGlen0(1)+zeros(MUAnew.Nnodes,1);
Fnew.n=F.n(1)+zeros(MUAnew.Nnodes,1);

[UserVar,Fnew]=GetGeometryAndDensities(UserVar,CtrlVar,MUAnew,Fnew,"-s-b-S-B-rho-");
%[UserVar,Fnew.s,Fnew.b,Fnew.S,Fnew.B,Fnew.rho,Fnew.rhow,Fnew.g]=DefineGeometryAndDensities(UserVar,CtrlVar,MUAnew,Fnew,"-s-b-S-B-rho-") ; 
Fnew.h=Fnew.s-Fnew.b;
Fnew.h0=Fnew.h; 

hf=Fnew.rhow.*(Fnew.S-Fnew.B)./Fnew.rho0 ;
Fnew.GF.node = HeavisideApprox(CtrlVar.kH,Fnew.h0-hf,CtrlVar.Hh0);  % 1 if grounded, 0 if afloat

%[UserVar,Fnew.C,Fnew.m,Fnew.q,Fnew.muk,Fnew.V0]=DefineSlipperyDistribution(UserVar,CtrlVar,MUAnew,Fnew);
% [UserVar,Fnew.C,Fnew.m]=DefineSlipperyDistribution(UserVar,CtrlVar,MUAnew,Fnew);
[UserVar,Fnew]=GetSlipperyDistribution(UserVar,CtrlVar,MUAnew,Fnew);

OutsideValues=[];
[RunInfo,Fnew.ub,Fnew.vb,Fnew.phi]=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,RunInfo,MUAold,MUAnew,OutsideValues,F.ub,F.vb,F.phi);

Fnew=PPFphi2F(CtrlVar,MUAnew,Fnew) ; 





end