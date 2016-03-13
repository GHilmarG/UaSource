function FEmeshCPT=CreateFEmeshCornerPointTriangulation(connectivity,coordinates)
    
    % FEmeshCPT=CreateFEmeshCornerPointTriangulation(connectivity,coordinates)
    %  creates a Triangulation representation
    % (see matlab documentation) from the corner nodes of the FE mesh
    %
    %
    
    
    [Nele,nod]=size(connectivity);
    
    if Nele==0
        FEmeshCPT=[];
        return
    end
    
    warning('off','MATLAB:triangulation:PtsNotInTriWarnId')
    
   
    
    switch nod
        case 3
            FEmeshCPT = triangulation(connectivity, coordinates(:,1),coordinates(:,2));
        case 6
            FEmeshCPT = triangulation(connectivity(:,1:2:5), coordinates(:,1),coordinates(:,2));
        case 10
            FEmeshCPT = triangulation(connectivity(:,1:3:7), coordinates(:,1),coordinates(:,2));
        otherwise
            error(' case not recognized')
    end
    warning('on','MATLAB:triangulation:PtsNotInTriWarnId')
    %
    %FEedges=edges(FEmeshCPT);
    %FEneighbors=neighbors(FEmeshCPT);
    %triplot(FEmeshCPT,'color','k') ; 
    %
end