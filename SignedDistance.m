function [DistSigned,isInside,isOnBoundary]=SignedDistance(Points,Boundary)

%%



DistSigned= pdist2(Boundary,Points,'euclidean','Smallest',1) ; 
DistSigned=DistSigned(:) ;

[isInside,isOnBoundary]=InsideOutside(Points,Boundary)  ; 
DistSigned(~isInside)=-DistSigned(~isInside);



return

end
