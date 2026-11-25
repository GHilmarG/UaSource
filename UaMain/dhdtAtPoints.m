function [hdot,Gu,Gv,Gh,Rhu,Rhv]=dhdtAtPoints(MUA,F,P,g)

%%
%
%   [hdot,Gu,Gv,Gh]=dhdtAtPoints(MUA,F,P)
%   [hdot,Gu,Gv,Gh,Rhu,Rhv]=dhdtAtPoints(MUA,F,P,g)
%
% Modelled rate of thickness change, and its derivatives with respect to the nodal values of u, v and h, at a set of points P,
% for example measurement locations (ShapeFunctionsAtPoints.m) or integration points (ShapeFunctionsAtIntegrationPoints.m):
%
% $$ \dot{h} = a - \frac{1}{\rho} \nabla \cdot (\rho h \mathbf{v}) = a - \left( \kappa_x u + h \partial_x u + \kappa_y v + h \partial_y v \right) ,
%    \qquad \kappa_x = \partial_x h + \frac{h}{\rho}\partial_x \rho , \quad \kappa_y = \partial_y h + \frac{h}{\rho}\partial_y \rho $$
%
% with a=F.as+F.ab (assumed independent of u, v and h), u=F.ub, v=F.vb, evaluated as in EvaluateJhdotAndDerivatives.m.
%
% Outputs (np = number of points):
%
%   hdot       : np x 1
%   Gu, Gv, Gh : np x Nnodes sparse Jacobians, d hdot_p / d u_j etc:
%
%     d hdot / d u_j = -( kappa_x phi_j + h d_x phi_j )
%     d hdot / d v_j = -( kappa_y phi_j + h d_y phi_j )
%     d hdot / d h_j = -( (d_x phi_j + phi_j d_x rho/rho) u + phi_j d_x u + (d_y phi_j + phi_j d_y rho/rho) v + phi_j d_y v )
%
%   Rhu, Rhv   : (only if g is given) Nnodes x Nnodes sparse, the g-weighted second derivatives
%
%     Rhu(j,i) = sum_p g_p d^2 hdot_p / (d h_j d u_i) ,  d^2 hdot/(d h_j d u_i) = -( (d_x phi_j + phi_j d_x rho/rho) phi_i + phi_j d_x phi_i )
%
%     and equivalently for v. These are the residual terms of the Hessian of a misfit with respect to (h,u) and (h,v), with g the
%     residuals weighted by the inverse error variances (and mapped back to the points), see dhdtObservationMisfit.m.
%
% hdot is linear in each of u, v and h separately, so there are no further second derivatives.
%
% see also: dhdtObservationMisfit.m, ShapeFunctionsAtPoints.m, ShapeFunctionsAtIntegrationPoints.m
%
%%

narginchk(3,4)

nod=MUA.nod ; N=MUA.Nnodes ; np=numel(P.EleID) ;
con=MUA.connectivity(P.EleID,:);      % np x nod

fun=P.fun ; Dx=P.dfdx ; Dy=P.dfdy ;

a=F.as+F.ab;
an=a(con) ; hn=F.h(con) ; un=F.ub(con) ; vn=F.vb(con) ; rn=F.rho(con) ;
if np==1 % with a single point, indexing returns a column, make sure these are rows
    an=an(:)' ; hn=hn(:)' ; un=un(:)' ; vn=vn(:)' ; rn=rn(:)' ;
end

ap=sum(fun.*an,2) ; h=sum(fun.*hn,2) ; u=sum(fun.*un,2) ; v=sum(fun.*vn,2) ; rho=sum(fun.*rn,2) ;
hx=sum(Dx.*hn,2) ; hy=sum(Dy.*hn,2) ; ux=sum(Dx.*un,2) ; vy=sum(Dy.*vn,2) ;
rx=sum(Dx.*rn,2) ; ry=sum(Dy.*rn,2) ;

kx=hx+h.*rx./rho ;
ky=hy+h.*ry./rho ;

hdot=ap-(kx.*u+h.*ux+ky.*v+h.*vy) ;

if nargout>1

    iRow=repmat((1:np)',1,nod);
    Px=Dx+fun.*(rx./rho) ;   % d_x phi_j + phi_j d_x rho/rho
    Py=Dy+fun.*(ry./rho) ;

    Gu=sparse(iRow(:),con(:),-(kx.*fun+h.*Dx),np,N);
    Gv=sparse(iRow(:),con(:),-(ky.*fun+h.*Dy),np,N);
    Gh=sparse(iRow(:),con(:),-(Px.*u+fun.*ux+Py.*v+fun.*vy),np,N);

end

if nargout>4

    if nargin<4 || isempty(g)
        error('dhdtAtPoints:gMissing','The weights g are required for the second-derivative terms Rhu and Rhv.')
    end
    g=g(:);

    Iu=zeros(np,nod,nod) ; Ju=Iu ; Xu=Iu ; Xv=Iu ;
    for j=1:nod          % thickness node
        for i=1:nod      % velocity node
            Iu(:,j,i)=con(:,j);
            Ju(:,j,i)=con(:,i);
            Xu(:,j,i)=-g.*(Px(:,j).*fun(:,i)+fun(:,j).*Dx(:,i));
            Xv(:,j,i)=-g.*(Py(:,j).*fun(:,i)+fun(:,j).*Dy(:,i));
        end
    end
    Rhu=sparse(Iu(:),Ju(:),Xu(:),N,N);
    Rhv=sparse(Iu(:),Ju(:),Xv(:),N,N);

end

end
