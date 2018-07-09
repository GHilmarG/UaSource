function [xc,yc]=FindCorners(xB,yB)

%
% A rough atempt at finding corners of a boundary
%
% Will only (possibly) work if there are quite a few boundary points between
% corners.
%


%K=convhull(x,y);
%K=boundary(x(:),y(:));
%plot(x(K),y(K),'b-*') ; axis equal ; 
%xB=x(K) ; yB=y(K);

xBdiff=gradient(xB) ; yBdiff=gradient(yB);

theta=atan2(xBdiff,yBdiff);
gradtheta=gradient(theta);

I=find(abs(gradtheta)>0.02);

plot(xB,yB,'x-c')
ind=I(2:3:end);
xc=xB(ind) ;  yc=yB(ind) ; 
hold on  ; plot(xc,yc,'o-r')

legend('Boundary','Corners')

end