
%%

Points=[0.5 0.5 ; 1 1 ; 2 2 ] ;
Boundary=[0 0 ; 1 0 ; 1 1 ; 0 1]; 

DistSigned=SignedDistance(Points,Boundary); 

DistSigned

%%
% Two boundaries with a NaNs between 
[X,Y]=meshgrid(-5:0.1:5,-5:0.1:5);
Points=[X(:), Y(:)] ; 
Boundary=[0 0 ; 1 0 ; 4 4 ; 3 3 ; -3 2 ; 0 1]; Boundary(end+1,:)=Boundary(1,:);
Npoints=100;  Boundary = interparc(Npoints,Boundary(:,1),Boundary(:,2),'linear');

Boundary2= [-4 -4 ; -2 -4 ; -2 -2 ; -4 -2 ; -4 -4] ; 
Npoints=100;  Boundary2 = interparc(Npoints,Boundary2(:,1),Boundary2(:,2),'linear');

Boundary=[ Boundary ; NaN NaN ; Boundary2] ; 


DistSigned=SignedDistance(Points,Boundary); 

figure
contourf(X,Y,reshape(DistSigned,size(X)),40) ; colorbar
hold on
hold on ; plot(Boundary(:,1),Boundary(:,2),'r')
axis equal