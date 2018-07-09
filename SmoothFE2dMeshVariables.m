
function [Zsmooth]=SmoothFE2dMeshVariables(coordinates,Zrough,Smoothness,dx,dy)
	
	x=coordinates(:,1) ; y=coordinates(:,2) ;
	
	xv=min(x):dx:max(x) ; yv=min(y):dy:max(y) ;
	
	[znew,xnew,ynew]=gridfit(x,y,Zrough,xv,yv,'smoothness',Smoothness);
	Zsmooth=interp2(xnew,ynew,znew,x,y);
	
end

