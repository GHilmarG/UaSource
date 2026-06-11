




function [xint,yint] = CalcIntegrationPointsCoordinates(MUA)


%%
%
% Calculates coordinates of element integration points.
%
% Calculates integration point values, and the (x,y) coordinates of the integration points
%
%
%   See also: Node2Int.m
%
%%



ndim=2;


coox=reshape(MUA.coordinates(MUA.connectivity,1),MUA.Nele,MUA.nod);
cooy=reshape(MUA.coordinates(MUA.connectivity,2),MUA.Nele,MUA.nod);


xint=zeros(MUA.Nele,MUA.nip) ; yint=zeros(MUA.Nele,MUA.nip);
for Iint=1:MUA.nip
    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    xint(:,Iint)=coox*fun;
    yint(:,Iint)=cooy*fun;

end


end