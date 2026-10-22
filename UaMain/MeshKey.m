function Key=MeshKey(MUA)

%% Key of a mesh, used to detect changes of the mesh
%
%   Key=MeshKey(MUA)
%
% (17 Oct 2026) Returns a structure with the fields
%
%   Key.Topology   [checksum of MUA.connectivity, size(MUA.connectivity), number of nodes]
%   Key.Geometry   [checksum of MUA.coordinates, size(MUA.coordinates)]
%
% computed with ArrayChecksum.m. Key.Topology changes with any change of the connectivity, including a renumbering of the nodes or the
% elements and a change of the order of the nodes within an element; Key.Geometry changes with any change of the coordinates. The keys
% are exact and platform independent, so they can be stored (eg in MUA.MeshKey, and in restart files) and compared later.
%
% CreateMUA.m and UpdateMUA.m store the key of the mesh in MUA.MeshKey, and UpdateMUA.m uses it to decide what has to be recalculated.
% AssemblyPatternCache.m stores Key.Topology in the assembly patterns (field TopologyKey).
%%

Key.Topology=[ArrayChecksum(MUA.connectivity) size(MUA.connectivity) size(MUA.coordinates,1)] ;
Key.Geometry=[ArrayChecksum(MUA.coordinates) size(MUA.coordinates)] ;

end
