function FEmeshTriRep=CreateFEmeshTriRep(connectivity,coordinates)

% [FEmeshTriRep]=CreateFEmeshTriRep(connectivity,coordinates)
%  creates a triangulation (see matlab documentation) from the corner nodes of the FE mesh
%
% Usefull for performing topological and geometric queries.
%
% See also TriFE

[Nele,nod]=size(connectivity);

if Nele==0
    FEmeshTriRep=[];
    return
end


warning('off','MATLAB:triangulation:PtsNotInTriWarnId')

switch nod
    case 3
        FEmeshTriRep = triangulation(connectivity, coordinates(:,1),coordinates(:,2));
    case 6
        FEmeshTriRep = triangulation(connectivity(:,1:2:5), coordinates(:,1),coordinates(:,2));
    case 10
        FEmeshTriRep = triangulation(connectivity(:,1:3:7), coordinates(:,1),coordinates(:,2));
    otherwise
        error(' case not recognized')
end

warning('on','MATLAB:triangulation:PtsNotInTriWarnId')
%warning('on','MATLAB:TriRep:PtsNotInTriWarnId')
% FEedges=edges(FEmeshTriRep);
% FEneighbors=neighbors(FEmeshTriRep)
% triplot(FEmeshTriRep,'color','k') ;



end