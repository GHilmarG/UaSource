




function [MeshDeriv,MeshDetJ]=CalcMuaMeshDerivatives(CtrlVar,MUA) 

% Note: CtrlVar is never used and can be left empty (2025 Dec)

narginchk(2,2)

dof=2; 
MeshDeriv=zeros(MUA.Nele,dof,MUA.nod,MUA.nip);
MeshDetJ=zeros(MUA.Nele,MUA.nip);



if MUA.Nele==0
        MeshDeriv=[];
        MeshDetJ=[];
        warning('CalcMeshDerivatives:EmptyMesh','Number of elements is zero!')
        return
end

% ndim=2;  points=sample('triangle',nip,ndim);

for Iint=1:MUA.nip
    [Deriv,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,MUA.points,Iint);
    MeshDeriv(:,:,:,Iint)=Deriv;
    MeshDetJ(:,Iint)=detJ;
end





end
