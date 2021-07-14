function [DTxy,tri,DTint,DTintTriInside,Xint,Yint,xint,yint,Iint]=TriangulationNodesIntegrationPoints(MUA)
    
    
    % Triangulation which is needed for plotting and interpolation purposes
    % DTxy           : Delaunay triangulation of all nodal points, not same as FE triangulation
    % tri            : 3-node triangulation
    % xint and yint  : integration point coordinates (Nele times nip)
    % Xint and Yint  : an array of unique integration points
    % DTint          : Delaunay triangulation of unique integration points (Xint,Yint)
    % DTintTriInside : triangulation only containing integration-point triangles with incenteres inside of MeshBoundaryCoordinates
    %
    % Iint           : index into integration-point variables needed for interpolation
    %
    % To interpolate the integration-point variable fint onto (x,y) use:
    %    temp=fint(:) ; temp=temp(Iint); fxy=Grid1toGrid2(DTint,temp,x,y)
    %
    
    
    
    x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
    % DTxy = DelaunayTri(x,y);
    DTxy = delaunayTriangulation(x,y);
    
    %ic=incenters(DTxy);
    %[cnNodes,on] = inpoly(ic,[x(BoundaryEdgeCornerNodes) y(BoundaryEdgeCornerNodes)]);
    %tri=DTxy.Triangulation(cnNodes,:);
    tri=TriFE(MUA.connectivity);
    
    %%
    
    
    
    if nargout > 2
        [xint,yint] = CalcIntegrationPointsCoordinates(MUA);
        
        % create vectors Xint and Yint of unique integration points and triangulise that set of points
        Xint=xint(:) ; Yint=yint(:); [~, Iint, ~] = unique([Xint Yint],'first','rows'); Iint = sort(Iint); Xint = Xint(Iint); Yint = Yint(Iint);
        %DTint = DelaunayTri(Xint,Yint);
        DTint = delaunayTriangulation(Xint,Yint);
        
        % get rid of triangles outside of the polygon define by MeshBoundaryCoordinates
        ic=incenters(DTint);
        [cnInt,on] = inpoly2(ic,[x(MUA.Boundary.EdgeCornerNodes) y(MUA.Boundary.EdgeCornerNodes)]);
        
        DTintTriInside=DTint.Triangulation(cnInt,:);
        
    end
    
end

