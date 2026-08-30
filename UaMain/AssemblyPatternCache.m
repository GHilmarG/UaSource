

function [uvAssemblyPattern,uvhAssemblyPattern]=AssemblyPatternCache(CtrlVar,MUA)

%%
% Precomputes, once per mesh, everything about the uv and/or the uvh
% assembly that depends only on MUA.connectivity.
%
%   [uvAssemblyPattern,uvhAssemblyPattern]=AssemblyPatternCache(CtrlVar,MUA)
%
% Which patterns are built is controlled by
%
%   CtrlVar.MUA.AssemblyPattern.uv  = true/false
%   CtrlVar.MUA.AssemblyPattern.uvh = true/false
%
% A scalar logical is also accepted for backwards compatibility and is
% taken to mean uv only.  Patterns that are not requested are returned
% empty.
%
% Usage in an assembly, for uv:
%
%   Xval = [d1d1(:) ; d2d2(:) ; d1d2(:) ; d2d1(:)];
%   vs   = accumarray(P.map,Xval,[P.nk 1]);
%   vs   = 0.5*(vs+vs(P.perm));                    % symmetrisation
%   Kuv  = sparse(P.i0,P.j0,vs,P.neq,P.neq);
%
% and for uvh, with the nine blocks in row-major order and no
% symmetrisation, since the uvh tangent is not symmetric:
%
%   Xval = [K11(:);K12(:);K13(:) ; K21(:);K22(:);K23(:) ; K31(:);K32(:);K33(:)];
%   vs   = accumarray(P.map,Xval,[P.nk 1]);
%   Kuvh = sparse(P.i0,P.j0,vs,P.neq,P.neq);
%
% Why this is faster
% ------------------
% The assembly hands sparse() d^2*nod^2*Nele unsorted triplets which
% collapse onto a much smaller set of nonzeros.  sparse() must sort them
% and sum the duplicates on every call, and the required permutation is
% identical every time because it depends only on the connectivity.  Doing
% that sort once turns the per-call work into one scatter-add plus a
% construction from already sorted, already unique triplets.
%
% What is shared between uv and uvh
% ---------------------------------
% Every block of both matrices has the same Nnodes x Nnodes sparsity, the
% nodal adjacency pattern of the mesh.  So the expensive step - the unique()
% over nod^2*Nele nodal keys - is done once here and both patterns are
% derived from it.  With nk1 nodal nonzeros, the uv matrix has 4*nk1 and
% the uvh matrix 9*nk1 nonzeros, exactly.
%
% The same nodal grouping underlies the level-set assemblies,
% MassMatrix2D1dof, StiffnessMatrix2D1dof and the dFuvd* family, so it is
% returned in the .Nodal field of each pattern in case those are converted
% later.
%
% Integer types
% -------------
% Linear keys are formed in uint64, never double: both are eight bytes and
% equally fast, but uint64 is exact to 2^64-1 whereas double is exact only
% to 2^53, and max(key)=neq^2.  That moves the point at which keys stop
% being exactly representable from about 47 million nodes to about 2.1
% billion.  MATLAB saturates rather than wraps on integer overflow, so
% there is an explicit guard.
%
% Stored index arrays are narrowed to uint32 where they fit and to uint64
% otherwise - never to double.  sparse has accepted integer subscripts
% since R2020a, and accumarray and ordinary indexing accept them too; the
% integer path is roughly twice as fast as the double path for all three.
%
% Element-array storage convention assumed is d**(:,Jnod,Inod), i.e.
% [Nele , j , i], so that d1d1(:) runs with the element index fastest,
% then Jnod, then Inod.
%
% Memory
% ------
% The dominant item is .map, of length d^2*nod^2*Nele.  For Nele=122850,
% nod=3 that is 40 MB for uv and 90 MB for uvh as uint32.  The shared
% .Nodal part is small (about 8 MB for the same mesh).
%
%%

narginchk(2,2)

[doUV,doUVH]=ParseRequest(CtrlVar);

uvAssemblyPattern=[];
uvhAssemblyPattern=[];

if ~doUV && ~doUVH
    return
end

Nele=MUA.Nele;
nod=MUA.nod;
Nnodes=MUA.Nnodes;

MaxNeq=floor(sqrt(double(intmax('uint64'))));
if 3*Nnodes > MaxNeq
    error('AssemblyPatternCache:tooLarge', ...
        'Nnodes=%d is too large for uint64 linear keys.',Nnodes)
end

conn=double(MUA.connectivity);

%% shared part: the nodal sparsity pattern

% row index for the (e,j,i) ordering: conn(e,Inod), broadcast over Jnod
I1=reshape(repmat(reshape(conn,Nele,1,nod),1,nod,1),[],1);
% column index: conn(e,Jnod), broadcast over Inod
J1=reshape(repmat(reshape(conn,Nele,nod,1),1,1,nod),[],1);

N64=uint64(Nnodes);
lin1=uint64(I1)+(uint64(J1)-1).*N64;

[u1,~,map1]=unique(lin1);

nk1=numel(u1);
i1=mod(u1-1,N64)+1;
% NB idivide with 'floor', not ./ : integer division in MATLAB rounds to
% nearest, which would be wrong here.
j1=idivide(u1-1,N64,'floor')+1;

Nodal.map=Narrow(map1,nk1);
Nodal.i1=Narrow(i1,Nnodes);
Nodal.j1=Narrow(j1,Nnodes);
Nodal.nk=nk1;
Nodal.Nnodes=Nnodes;
Nodal.Nele=Nele;
Nodal.nod=nod;

%% derive the requested patterns

if doUV
    % block order must match Xval=[d1d1(:) ; d2d2(:) ; d1d2(:) ; d2d1(:)]
    blocks=[1 1 ; 2 2 ; 1 2 ; 2 1];
    uvAssemblyPattern=BuildPattern(conn,i1,j1,map1,nk1,Nnodes,Nele,nod,2,blocks,true);
    uvAssemblyPattern.Nodal=Nodal;
end

if doUVH
    % nine blocks in row-major order, matching the loop in uvhMatrixAssembly
    blocks=[1 1 ; 1 2 ; 1 3 ; 2 1 ; 2 2 ; 2 3 ; 3 1 ; 3 2 ; 3 3];
    uvhAssemblyPattern=BuildPattern(conn,i1,j1,map1,nk1,Nnodes,Nele,nod,3,blocks,false);
    uvhAssemblyPattern.Nodal=Nodal;
end

end




function P=BuildPattern(conn,i1,j1,map1,nk1,Nnodes,Nele,nod,d,blocks,wantPerm)

% Builds the pattern for a d-field system whose blocks, in the order the
% assembly concatenates them, are listed in blocks as [rowfield colfield].

nb=size(blocks,1);
neq=d*Nnodes;
neq64=uint64(neq);
n1t=nod*nod*Nele;      % triplets per block
nk=nb*nk1;             % blocks are disjoint, so no duplicates across them

% global column-major key of every nonzero, block by block
g=zeros(nk,1,'uint64');
for c=1:nb
    a=uint64((blocks(c,1)-1)*Nnodes);
    b=uint64((blocks(c,2)-1)*Nnodes);
    g((c-1)*nk1+1:c*nk1)=(i1+a)+((j1+b)-1).*neq64;
end

[gs,ord]=sort(g);

% pos(k) is where entry k of the block-stacked list ends up once sorted
pos=zeros(nk,1,'uint32');
if nk>intmax('uint32')
    pos=zeros(nk,1,'uint64'); pos(ord)=uint64(1:nk);
else
    pos(ord)=uint32(1:nk);
end

i0=mod(gs-1,neq64)+1;
j0=idivide(gs-1,neq64,'floor')+1;

% map from every triplet to its position in the sorted nonzero list
if nk>intmax('uint32')
    map=zeros(nb*n1t,1,'uint64');
else
    map=zeros(nb*n1t,1,'uint32');
end
for c=1:nb
    posc=pos((c-1)*nk1+1:c*nk1);
    map((c-1)*n1t+1:c*n1t)=posc(map1);
end

P.map=map;
P.i0=Narrow(i0,neq);
P.j0=Narrow(j0,neq);
P.nk=nk;
P.neq=neq;
P.blocks=blocks;
P.d=d;

% transpose permutation, for folding (K+K')/2 into the value array.  Only
% meaningful where the assembled matrix is symmetric, i.e. uv.
if wantPerm
    linT=j0+(i0-1).*neq64;
    [tf,perm]=ismember(linT,gs);
    if ~all(tf)
        error('AssemblyPatternCache:notSymmetric', ...
            'The assembly sparsity pattern is not structurally symmetric.')
    end
    P.perm=Narrow(perm,nk);
else
    P.perm=[];
end

% indices for the residual / external-force assembly: Inod outer, field inner
iR=zeros(nod*Nele*d,1);
istak=0;
for Inod=1:nod
    for a=1:d
        iR(istak+1:istak+Nele)=conn(:,Inod)+(a-1)*Nnodes;
        istak=istak+Nele;
    end
end
P.iR=Narrow(iR,neq);

P.Nele=Nele;
P.nod=nod;
P.Nnodes=Nnodes;

end




function v=Narrow(v,bound)

% uint32 where it fits, uint64 otherwise - never double.

if bound<=intmax('uint32')
    v=uint32(v);
else
    v=uint64(v);
end

end




function [doUV,doUVH]=ParseRequest(CtrlVar)

doUV=false; doUVH=false;

if ~isfield(CtrlVar,'MUA') || ~isfield(CtrlVar.MUA,'AssemblyPattern')
    return
end

A=CtrlVar.MUA.AssemblyPattern;

if isstruct(A)
    if isfield(A,'uv')  ; doUV =logical(A.uv)  ; end
    if isfield(A,'uvh') ; doUVH=logical(A.uvh) ; end
else
    % backwards compatibility: a scalar logical means uv only
    doUV=logical(A);
end

end
