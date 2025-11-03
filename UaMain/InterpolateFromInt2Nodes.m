function fnod=InterpolateFromInt2Nodes(fint,coordinates,xint,yint,DTint,Iint)

% now better to use InterpolateIntVariables.m
    
if nargin==4
    
    
    Xint=xint(:) ; Yint=yint(:);
    [~, Iint, ~] = unique([Xint Yint],'first','rows');
    Iint = sort(Iint); Xint = Xint(Iint); Yint = Yint(Iint);
    DTint = DelaunayTri(Xint,Yint);
    %ic=incenters(DTint); [cn,on] = inpoly(ic,MeshBoundaryCoordinates); DTintTriInside=DTint.Triangulation(cn,:);
    
end

ftemp=fint(:) ; ftemp=ftemp(Iint);

F = TriScatteredInterp(DTint,real(ftemp),'natural');
fnod = F(coordinates(:,1),coordinates(:,2));

% there will be some nodes along the boundary that are not within the convex set
% ind=find(isnan(w))

if any(isnan(fnod))
    
    
    ind=find(isnan(fnod));
    
    nn=nearestNeighbor(DTint,coordinates(ind,1),coordinates(ind,2));
    fnod(ind)=ftemp(nn);
    
end


end

