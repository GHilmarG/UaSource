function P=ShapeFunctionsAtIntegrationPoints(CtrlVar,MUA)

%%
%
%   P=ShapeFunctionsAtIntegrationPoints(CtrlVar,MUA)
%
% Form functions, their x and y derivatives, and the quadrature weights at all integration points of the mesh, in the same format as
% returned by ShapeFunctionsAtPoints.m. The integration points are ordered element by element within each integration point, i.e. the
% point (element e, integration point q) has the index e+(q-1)*MUA.Nele.
%
%   P.EleID   : (Nele*nip) x 1 element number
%   P.fun     : (Nele*nip) x nod form functions
%   P.dfdx    : (Nele*nip) x nod x-derivatives of the form functions
%   P.dfdy    : (Nele*nip) x nod y-derivatives of the form functions
%   P.w       : (Nele*nip) x 1 quadrature weights (detJ times the weight of the integration rule), sum(P.w) is the area of the mesh
%   P.x, P.y  : (Nele*nip) x 1 coordinates of the integration points
%
% see also: ShapeFunctionsAtPoints.m, dhdtAtPoints.m
%
%%

narginchk(2,2)

if ~isfield(MUA,'Deriv') || isempty(MUA.Deriv) || ~isfield(MUA,'DetJ') || isempty(MUA.DetJ)
    [MUA.Deriv,MUA.DetJ]=CalcMuaMeshDerivatives(CtrlVar,MUA);
end

Nele=MUA.Nele ; nod=MUA.nod ; nip=MUA.nip ; np=Nele*nip ;

P.EleID=repmat((1:Nele)',nip,1);
P.fun=zeros(np,nod) ; P.dfdx=zeros(np,nod) ; P.dfdy=zeros(np,nod) ;
P.w=zeros(np,1) ; P.x=zeros(np,1) ; P.y=zeros(np,1) ;

xnod=reshape(MUA.coordinates(MUA.connectivity,1),Nele,nod);
ynod=reshape(MUA.coordinates(MUA.connectivity,2),Nele,nod);

for Iint=1:nip
    I=(Iint-1)*Nele+(1:Nele);
    fun=shape_fun(Iint,2,nod,MUA.points);       % nod x 1
    P.fun(I,:)=repmat(fun',Nele,1);
    P.dfdx(I,:)=reshape(MUA.Deriv(:,1,:,Iint),Nele,nod);
    P.dfdy(I,:)=reshape(MUA.Deriv(:,2,:,Iint),Nele,nod);
    P.w(I)=MUA.DetJ(:,Iint)*MUA.weights(Iint);
    P.x(I)=xnod*fun;
    P.y(I)=ynod*fun;
end

end
