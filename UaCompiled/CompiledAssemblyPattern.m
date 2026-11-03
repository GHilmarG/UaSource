function C=CompiledAssemblyPattern(P)

%% Per-mesh data for the compiled assembly (uv and uvh), derived from the cached assembly pattern P (MUA.uvAssemblyPattern or MUA.uvhAssemblyPattern)
%
%   invPtr, invIdx : for each nonzero k of the sparse matrix (in CSC order), the positions in the element-major array Xmat(L,e) of the
%                    element-matrix values that are summed into it, in the same order as accumarray(P.map,...) sums them
%   perm           : for the symmetrisation of the nonzeros (uv: P.perm; uvh: empty)
%   rPtr, rIdx     : for each degree of freedom r (row of the residual), the positions in the element-major array of the element
%                    residuals Rmat(L,e), L = i + nod*(c-1) (local node i, component c), that are summed into it, in the same order as
%                    accumarray(P.iR,...) sums them (31 Oct 2026: the residual is summed in the kernel)
%   jc, ir         : zero-based column pointers and row indices of the matrix (SparseFromCSC)
%   diagIdx        : positions of the diagonal entries among the nonzeros (for the NaN check)
%
% The number of blocks of the element matrix (4 for uv, 9 for uvh) follows from the size of P.map. Kept in a persistent cache with one
% entry per assembly (P.d=2: uv, P.d=3: uvh), keyed by the topology key of the mesh, so it is only recomputed when the mesh changes, and
% uv and uvh assemblies do not evict each other.
%
% (30 Oct 2026, generalised from uvhCSCPattern.m)
%%

persistent Cache
if isempty(Cache) ; Cache=cell(1,3) ; end
slot=double(P.d);
if ~isempty(Cache{slot}) && isequal(Cache{slot}.Key,P.TopologyKey) && Cache{slot}.nk==P.nk && Cache{slot}.Nele==P.Nele
    C=Cache{slot}.C; return
end

Nele=double(P.Nele); nod=double(P.nod);
nb=numel(P.map)/(Nele*nod*nod);                        % number of blocks of the element matrix: 4 (uv) or 9 (uvh)
nL=nb*nod*nod;
[~,t]=sort(P.map);                                     % stable: positions t in Xval (block-major), grouped by nonzero, increasing t
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
% residual: P.iR gives the row of each entry of Tval=reshape(permute(cat(3,Tx,Ty[,Th]),[1 3 2]),[],1), ie position
% t = e + Nele*(c-1) + Nele*d*(i-1) (element e, component c, local node i)
d=double(P.d);
[~,tr]=sort(P.iR);                                     % stable: positions in Tval, grouped by row, increasing t
tr=double(tr);
er=mod(tr-1,Nele)+1;
qr=floor((tr-1)/Nele);
cr=mod(qr,d)+1; ir=floor(qr/d)+1;
C.rIdx=uint32(ir+nod*(cr-1)+d*nod*(er-1));             % position in the element-major array Rmat(L,e)
C.rPtr=uint32([0;cumsum(accumarray(double(P.iR),1,[P.neq 1]))]);
if isfield(P,"perm") && ~isempty(P.perm)
    C.perm=uint32(P.perm);
else
    C.perm=zeros(0,1,"uint32");
end
Cache{slot}=struct("C",C,"Key",P.TopologyKey,"nk",P.nk,"Nele",P.Nele);
end
