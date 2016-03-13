function [coordinates,connectivity]=UserMeshModifications(CtrlVar,coordinates,connectivity)

% This m-file is called directly after each meshing step to allow for 
% further modifications of the FE mesh
%
% [coordinates,connectivity]=UserMeshModifications(CtrlVar,coordinates,connectivity)
%
% Example: Delete selected elements from the mesh
%     
%    DeactivatedElements=[10 ; 11] ;  % list of elements to be deactivated
%    [coordinates,connectivity]=DeactivateElements(CtrlVar,iDeactivatedElements,coordinates,connectivity);
%
%
%
%

end 