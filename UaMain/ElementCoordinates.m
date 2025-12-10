function [xEle,yEle,DTele]=ElementCoordinates(connectivity,coordinates)

% calculates centre points of triangles, and the triangularisation of
% those centre points


[Nele,nod]=size(connectivity);



xEle=mean(reshape(coordinates(connectivity,1),Nele,nod),2);
yEle=mean(reshape(coordinates(connectivity,2),Nele,nod),2);

if nargout>2
    DTele= delaunayTriangulation(xEle,yEle);
    %DTele = DelaunayTri(xEle,yEle);
end
%TRIele=DTele.Triangulation;

end



