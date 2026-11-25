function P=ShapeFunctionsAtPoints(CtrlVar,MUA,x,y)

%%
%
%   P=ShapeFunctionsAtPoints(CtrlVar,MUA,x,y)
%
% Form functions, and their derivatives with respect to x and y, at arbitrary locations (x,y) within the finite-element mesh.
%
% On return, P is a structure with the fields
%
%   P.Inside  : numel(x) x 1 logical, true if (x(i),y(i)) lies within the mesh. All other fields refer to the inside locations only.
%   P.EleID   : np x 1 element containing each (inside) location
%   P.fun     : np x MUA.nod values of the form functions
%   P.dfdx    : np x MUA.nod derivatives of the form functions with respect to x
%   P.dfdy    : np x MUA.nod derivatives of the form functions with respect to y
%   P.x, P.y  : np x 1 the (inside) locations
%
% so that, for a nodal field f, the value and the x-derivative at the inside locations are
%
%   sum(P.fun.*f(MUA.connectivity(P.EleID,:)),2)   and   sum(P.dfdx.*f(MUA.connectivity(P.EleID,:)),2)
%
% The enclosing elements are found using BuildNode2DataMap.m, and the form functions evaluated using sr_shape_fun.m at the barycentric
% coordinates (see BuildNode2DataMap.m). All elements are straight sided, and the barycentric coordinates are therefore affine functions
% of x and y, with constant gradients within each element. The derivatives of the form functions are calculated exactly (to round-off)
% by complex-step differentiation of sr_shape_fun.m along these gradients.
%
% see also: BuildNode2DataMap.m, ShapeFunctionsAtIntegrationPoints.m, dhdtAtPoints.m
%
%%

narginchk(4,4)

x=x(:) ; y=y(:) ;

[~,Inside,ID,MUA]=BuildNode2DataMap(CtrlVar,MUA,x,y);

P.Inside=Inside;
P.x=x(Inside);
P.y=y(Inside);
P.EleID=ID(Inside);

TR=MUA.TR;
B=cartesianToBarycentric(TR,P.EleID,[P.x P.y]);

V=TR.ConnectivityList(P.EleID,:);
x1=TR.Points(V(:,1),1) ; y1=TR.Points(V(:,1),2) ;
x2=TR.Points(V(:,2),1) ; y2=TR.Points(V(:,2),2) ;
x3=TR.Points(V(:,3),1) ; y3=TR.Points(V(:,3),2) ;

TwoA=(x2-x1).*(y3-y1)-(x3-x1).*(y2-y1);         % twice the (signed) element area

dBdx=[y2-y3 , y3-y1 , y1-y2]./TwoA;              % gradients of the barycentric coordinates
dBdy=[x3-x2 , x1-x3 , x2-x1]./TwoA;

P.fun=sr_shape_fun(B,MUA.nod);

hcs=1e-20;   % complex step
P.dfdx=imag(sr_shape_fun(B+1i*hcs*dBdx,MUA.nod))/hcs;
P.dfdy=imag(sr_shape_fun(B+1i*hcs*dBdy,MUA.nod))/hcs;

end
