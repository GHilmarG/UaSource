

function [Qn,Qt,qn,qt,xc,yc,normal,t]=PathIntegral(CtrlVar,x,y,qx,qy)

%%
%
% Simple utility to calculate normal and tangential flux for each line segment defined by the (x,y) array.
%
% Also returns the total summed values
%
%
% Qn    : total summed flux normal to the curve
% Qt    : total summed flux tangential to the curve
%
% qn    : normal flux for each line segment
% qt    : tangential flux for each line segment
%
%
% qx and qy can be either:
%         1) vectors with the q values at the locations (x,y)
%         2) interpolants providing qx and qy values at (x,y)
%
%
%
%
% Example:
%
% 
% BoundaryNodes=MUA.Boundary.EdgeCornerNodes ;      % Here the boundary is defined as the boundary of the mesh
%                                                   % But note that this only includes corner nodes
% BoundaryNodes=[BoundaryNodes;BoundaryNodes(1)] ;  % Make sure to close the loop
% BoundaryNodes=flipud(BoundaryNodes) ;             % Depending on the desired orientation of the normals, it might be
%                                                   % required to flip the orientation of the curve.
% 
% [Qn,Qt,qn,qt,xc,yc,normal]=PathIntegral(CtrlVar,F.x(BoundaryNodes),F.y(BoundaryNodes),qxN(BoundaryNodes),qyN(BoundaryNodes));
% 
% Example using interpolants
%
% Fqx=scatteredInterpolant(F.x,F.y,qx) ;  Fqy=scatteredInterpolant(F.x,F.y,qy) ;
% [Qn,Qt,qn,qt,xc,yc,normal]=PathIntegral(CtrlVar,F.x(BoundaryNodes),F.y(BoundaryNodes),Fqx,Fqy);
%
%
% Note: This simple routine just takes the (x,y) values provided, and does not subdivide the curve further. To create set of
% points at fixed distance along a cure, consider uisng "interparc.m" by John D'Errico which is included in the UaSource
% directory, or can be donwloaded from Matlab central file exchange.
%
% Npoints=300 ; pt=interparc(Npoints,x,y,'linear');  figure(1) ; plot(pt(:,1),pt(:,2),"o-") ; axis equal
%
%




isInterpolant=isa(qx,"scatteredInterpolant") || isa(qx,"griddedInterpolant") ; 


N=numel(x);

qn=nan(N-1,1);
qt=nan(N-1,1);
xc=nan(N-1,1);
yc=nan(N-1,1);
normal=nan(N-1,2); 

for I=1:N-1

    dx=x(I+1)-x(I);
    dy=y(I+1)-y(I);
  
    xc(I)=x(I)+dx/2;  
    yc(I)=y(I)+dy/2;

    ds=sqrt(dx.*dx+dy.*dy);

    n=[dy ; -dx]./ds;
    t=[dx ; dy]./ds ;

    if isInterpolant
        qc=[qx(xc(I),yc(I)) ; qy(xc(I),yc(I))] ;
    else
        qc=[(qx(I+1)+qx(I))/2 ; (qy(I+1)+qy(I))/2] ;
    end

    % q=ds*[(qx(I+1)+qx(I))/2 ; (qy(I+1)+qy(I))/2] ;
    q=ds*qc;

    qn(I)=q'*n;  % Normal flux
    qt(I)=q'*t ;  % Tangential flux

    normal(I,:)=n'; 

end

Qn=sum(qn);
Qt=sum(qt);


end