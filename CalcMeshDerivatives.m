function [MeshDeriv,MeshDetJ]=CalcMeshDerivatives(CtrlVar,connectivity,coordinates)


dof=2; [Nele,nod]=size(connectivity);
nip=NrOfIntegrationPoints(CtrlVar);
MeshDeriv=zeros(Nele,dof,nod,nip);
MeshDetJ=zeros(Nele,nip);

if CtrlVar.InfoLevelCPU>=10 ;   tDeriv=tic; end
if CtrlVar.InfoLevel>=10;
    fprintf('CalcMeshDerivatives: calculating mesh derivatives \n ')
end

parfor Iint=1:nip
    [Deriv,detJ]=derivVector(coordinates,connectivity,nip,Iint);
    MeshDeriv(:,:,:,Iint)=Deriv;
    MeshDetJ(:,Iint)=detJ;
end


if CtrlVar.InfoLevelCPU>=10 ;  fprintf(CtrlVar.fidlog,'CalcMeshDerivatives in %-g sec. \n',toc(tDeriv)) ; end


end
