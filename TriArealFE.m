function A=TriAreaFE(coordinates,connectivity)
    
    % calculates the area of triangles in a FE mesh given coordinates and connectivity
	% assumes straigth sides
	%
    [Nele,nod]=size(connectivity);
    xnod=reshape(coordinates(connectivity,1),Nele,nod);
    ynod=reshape(coordinates(connectivity,2),Nele,nod);
    
    switch nod
        case 3
            A=TriArea(xnod,ynod);
        case 6
            A=TriArea(xnod(:,[1:2:6]),ynod(:,[1:2:6]));
        case 10
            A=TriArea(xnod(:,[1 4 7]),ynod(:,[1 4 7]));
        otherwise
            fprintf(' case not implemented ')
    end
    
    
end
