function C=uvhCSCPattern(P)
% Per-mesh data for the compiled uvh assembly, derived from the cached assembly pattern P=MUA.uvhAssemblyPattern:
%   invPtr, invIdx : for each nonzero k of the sparse matrix (in CSC order), the positions in the element-major array Xmat(L,e) of the
%                    element-matrix values that are summed into it (GatherSumCSC), in the same order as accumarray(P.map,...) sums them
%   jc, ir         : zero-based column pointers and row indices of the matrix (SparseFromCSC)
%   diagIdx        : positions of the diagonal entries among the nonzeros (for the NaN check)
% Kept in a persistent cache, keyed by the topology key of the mesh, so it is only recomputed when the mesh changes.
persistent Cache
if ~isempty(Cache) && isequal(Cache.Key,P.TopologyKey) && Cache.nk==P.nk && Cache.Nele==P.Nele
    C=Cache.C; return
end
Nele=double(P.Nele); nod=double(P.nod); nL=9*nod*nod;
[~,t]=sort(P.map);                                     % stable: positions t in Xval=[Kxu(:);...;Khh(:)] (block-major), grouped by nonzero, increasing t
t=double(t);
e=mod(t-1,Nele)+1;                                     % element
L=floor((t-1)/Nele)+1;                                 % position within the element matrix (block, j, i)
C.invIdx=uint32(L+nL*(e-1));                           % position in the element-major array Xmat(L,e)
cnt=accumarray(double(P.map),1,[P.nk 1]);
C.invPtr=uint32([0;cumsum(cnt)]);
cj=accumarray(double(P.j0),1,[P.neq 1]);
C.jc=uint64([0;cumsum(cj)]);
C.ir=uint64(P.i0)-1;
C.diagIdx=uint32(find(P.i0==P.j0));
C.neq=P.neq;
Cache.C=C; Cache.Key=P.TopologyKey; Cache.nk=P.nk; Cache.Nele=P.Nele;
end
