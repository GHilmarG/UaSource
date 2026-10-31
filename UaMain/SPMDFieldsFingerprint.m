function Key=SPMDFieldsFingerprint(MUA,nW,varargin)

%% Fingerprint of the inputs that the parallel (spmd) assemblies keep on the workers between calls
%
%   Key=SPMDFieldsFingerprint(MUA,nW,F,Exclude)
%   Key=SPMDFieldsFingerprint(MUA,nW,F0,Exclude0,F1,Exclude1,...)
%
% (17 Oct 2026) Used by uvMatrixAssemblySSTREAM_SPMD.m and uvhMatrixAssemblySSTREAM_SPMD.m to decide whether the fields kept on the
% workers have to be sent again. Key consists of the numbers of nodes and elements, the key of the mesh (MUA.MeshKey, see MeshKey.m), the
% number of workers nW, and, for each pair (F,Exclude), weighted sums of all numerical fields of F except those listed in Exclude (the
% fields that change between calls and are sent in each call). Fields that are structures (eg GF, LSFMask) contribute weighted sums over
% their numerical fields, and for each field the number of non-finite values is included as well.
%
% The weighted sums are fast, and a change of a field that leaves both its weighted sum and its number of non-finite values unchanged is
% extremely unlikely. In any case, the worst case of a detected change is an unnecessary resend, not stale data on the workers.
%
% (Previously the local functions StaticFingerprint, FieldSums and WeightedSum of uvhMatrixAssemblySSTREAM_SPMD.m; the key of the uvh
% assembly is unchanged.)
%
%%

if isfield(MUA,"MeshKey") && isstruct(MUA.MeshKey) && isfield(MUA.MeshKey,"Topology") && isfield(MUA.MeshKey,"Geometry")
    MK=[MUA.MeshKey.Topology MUA.MeshKey.Geometry] ;
else
    MK=MeshKey(MUA) ; MK=[MK.Topology MK.Geometry] ;
end

Key=[MUA.Nnodes MUA.Nele MK nW] ;
for k=1:2:numel(varargin)
    Key=[Key FieldSums(varargin{k},string(varargin{k+1}))] ; %#ok<AGROW>
end

end

%%
function s=FieldSums(F,Exclude)
p=setdiff(string(properties(F)),Exclude,"stable") ;
s=zeros(1,2*numel(p)) ;
for k=1:numel(p)
    x=F.(p(k)) ;
    if isstruct(x)          % eg LSFMask, GF: weighted sum over its numerical fields
        fn=fieldnames(x) ;
        for i=1:numel(fn)
            y=x.(fn{i}) ;
            if (isnumeric(y) || islogical(y)) && ~isempty(y)
                [a,b]=WeightedSum(y) ; s(2*k-1)=s(2*k-1)+i*a ; s(2*k)=s(2*k)+b ;
            end
        end
    elseif (isnumeric(x) || islogical(x)) && ~isempty(x)
        [s(2*k-1),s(2*k)]=WeightedSum(x) ;
    end
end
end

%%
function [v,n]=WeightedSum(x)
% weighted sum of the finite values, and the number of non-finite values
%
% (31 Oct 2026) The weights are kept between calls (persistent, by length), and the sum is a dot product. If it is finite, x has
% no non-finite values (the weights are >= 1); only otherwise are the non-finite values counted and left out. About 15 ms faster per
% fingerprint of a UaFields object (Greenland, 56,000 nodes). The values of the fingerprint differ from those of the previous version,
% which does not matter, as they are only compared with each other within a session.
persistent W
x=x(:) ; if ~isa(x,"double") || issparse(x) ; x=double(full(x)) ; end
N=numel(x) ;
if numel(W)<N
    W=1+mod((1:max(N,2*numel(W)))',97)/97 ;
end
v=W(1:N).'*x ; n=0 ;
if ~isfinite(v)
    nf=~isfinite(x) ; n=nnz(nf) ; x(nf)=0 ;
    v=W(1:N).'*x ;
end
end
