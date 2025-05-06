



function [MeshDeriv,MeshDetJ]=CalcMeshDerivatives(connectivity,coordinates,nip,points)

% Note: Do not use. 
% 
% Use 
% 
%  CalcMuaMeshDerivatives.m
%
% instead.

error("not longer used. Use CalcMuaMeshDerivatives.m instead.")

narginchk(4,4)

dof=2; [Nele,nod]=size(connectivity);
MeshDeriv=zeros(Nele,dof,nod,nip);
MeshDetJ=zeros(Nele,nip);



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




end
