function [UserVar,coordinates,connectivity]=GetMeshModifications(UserVar,CtrlVar,coordinates,connectivity)

N=nargout('DefineMeshModifications');

switch N
    
    case 2
        
        [coordinates,connectivity]=DefineMeshModifications(CtrlVar,coordinates,connectivity);
        
    case 3

        [UserVar,coordinates,connectivity]=DefineMeshModifications(UserVar,CtrlVar,coordinates,connectivity);
        
end

end
