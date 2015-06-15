function FEmeshTriRep=CreateFEmeshTriRep(connectivity,coordinates)
    
    % [FEmeshTriRep]=CreateFEmeshTriRep(connectivity,coordinates)
    %  creates a TriRep (see matlab documentation) from the corner nodes of the FE mesh
    %
    %
    [Nele,nod]=size(connectivity);
    
    if Nele==0
        FEmeshTriRep=[];
        return
    end
    
    warning('off','MATLAB:TriRep:PtsNotInTriWarnId')
    
    
    
    switch nod
        case 3
            FEmeshTriRep = TriRep(connectivity, coordinates(:,1),coordinates(:,2));
        case 6
            FEmeshTriRep = TriRep(connectivity(:,1:2:5), coordinates(:,1),coordinates(:,2));
        case 10
            FEmeshTriRep = TriRep(connectivity(:,1:3:7), coordinates(:,1),coordinates(:,2));
        otherwise
            error(' case not recognized')
    end
    warning('on','MATLAB:TriRep:PtsNotInTriWarnId')
    % FEedges=edges(FEmeshTriRep);
    % FEneighbors=neighbors(FEmeshTriRep)
    % triplot(FEmeshTriRep,'color','k') ; 
end