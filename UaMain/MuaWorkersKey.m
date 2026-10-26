function Key=MuaWorkersKey(MeshKey,MUA,nWorkers)

%% Key of the worker data (MUA.workers), used by UpdateMUA.m to decide if the worker data have to be built (anew)
%
%   Key=MuaWorkersKey(MeshKey,MUA,nWorkers)
%
% The worker data are the parts of the mesh, held by each worker of the parallel pool for the spmd assembly (see BuildMuaWorkers.m).
% They depend on the mesh (MeshKey, see MeshKey.m), the element type and the quadrature rule, and the number of workers. These are
% recorded in MUA.workersKey when the worker data are built (CreateMUA.m, UpdateMUA.m), and compared with the current ones by
% UpdateMUA.m.
%
% The contents of a Composite can only be checked by communicating with the workers, so this key is kept on the client instead.
% Whether the Composite still belongs to the current parallel pool is checked separately, with all(exist(MUA.workers)): this is false
% for all workers if the pool has been deleted, or replaced by a new one (even of the same size).
%
% See also: UpdateMUA, BuildMuaWorkers, MeshKey
%
% (26 Oct 2026)
%%

Key=struct("Mesh",MeshKey,"nod",MUA.nod,"nip",MUA.nip,"points",MUA.points,"weights",MUA.weights,"nWorkers",nWorkers);

end
