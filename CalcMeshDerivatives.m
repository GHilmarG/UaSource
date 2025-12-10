function [MeshDeriv,MeshDetJ]=CalcMeshDerivatives(CtrlVar,connectivity,coordinates,nip,points)


narginchk(5,5)

dof=2; [Nele,nod]=size(connectivity);
MeshDeriv=zeros(Nele,dof,nod,nip);
MeshDetJ=zeros(Nele,nip);

if CtrlVar.InfoLevelCPU>=10 ;   tDeriv=tic; end
if CtrlVar.InfoLevel>=10
    fprintf('CalcMeshDerivatives: calculating mesh derivatives \n ')
end

if Nele==0
        MeshDeriv=[];
        MeshDetJ=[];
        warning('CalcMeshDerivatives:EmptyMesh','Number of elements is zero!')
        return
end

% ndim=2;  points=sample('triangle',nip,ndim);

for Iint=1:nip
    [Deriv,detJ]=derivVector(coordinates,connectivity,nip,points,Iint);
    MeshDeriv(:,:,:,Iint)=Deriv;
    MeshDetJ(:,Iint)=detJ;
end


if CtrlVar.InfoLevelCPU>=10 ;  fprintf(CtrlVar.fidlog,'CalcMeshDerivatives in %-g sec. \n',toc(tDeriv)) ; end


end
