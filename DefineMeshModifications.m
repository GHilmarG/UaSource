function [coordinates,connectivity]=DefineMeshModifications(CtrlVar,coordinates,connectivity)

%%
% This m-file is called directly after each meshing step to allow for 
% further modifications of the FE mesh.
%
% [coordinates,connectivity]=UserMeshModifications(CtrlVar,coordinates,connectivity)
%
% A typical use is to deactivate some elements. 
% When deactivating elements, use the m-file 
%   DeactivateElements(CtrlVar,iDeactivatedElements,coordinates,connectivity);
% from within this m-file.
%
% Example: Deactivate selected elements from the mesh
%     
%    DeactivatedElements=[10 ; 11] ;  % list of elements to be deactivated
%    [coordinates,connectivity]=DeactivateElements(CtrlVar,iDeactivatedElements,coordinates,connectivity);
%
% Example: Deactivate elements whose centrepoints are within a given polygon: 
%     
%     load('NameOfFileWithPolygonEdges',xPoly,yPoly);
%     xEle=Nodes2EleMean(connectivity,coordinates(:,1));
%     yEle=Nodes2EleMean(connectivity,coordinates(:,2));
%     [iDeactivatedElements,on] = inpoly([xEle yEle],[xPoly yPoly]);
%     [coordinates,connectivity]=DeactivateElements(CtrlVar,iDeactivatedElements,coordinates,connectivity);
%
% See also DeactivateElements, inpoly
%%

end 