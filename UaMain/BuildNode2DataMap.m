


function [O,Inside,ID,MUA]=BuildNode2DataMap(CtrlVar,MUA,x,y)

%%
%
%   [O,Inside,ID,MUA]=BuildNode2DataMap(CtrlVar,MUA,x,y)
%
% Builds the linear mapping (observation) operator, O, that maps nodal values of the
% two-dimensional finite-element mesh MUA onto an arbitrary set of (x,y) locations,
% using the FE form (shape) functions.
%
% If F is any nodal field (a MUA.Nnodes x 1 vector), then
%
%   Fmeas = O * F
%
% is a m x 1 vector containing the FE estimate of that field at the m locations (x,y).
%
%
%% Inputs
%
%   CtrlVar     : Ua control variable structure. Only CtrlVar.InfoLevel is used, and the
%                 call also works with CtrlVar=[].
%
%   MUA         : Ua mesh structure. Uses MUA.connectivity, MUA.coordinates, MUA.Nnodes,
%                 MUA.Nele, MUA.nod and, if available, MUA.TR and MUA.EleAreas.
%
%   x , y       : vectors with the x and y locations onto which the nodal values are to
%                 be mapped, for example the locations of measurements. These are
%                 arbitrary locations and need not have anything to do with the nodal
%                 locations of the FE mesh. Both are reshaped internally into column
%                 vectors, so any orientation is accepted.
%
%
%% Outputs
%
%   O           : m x MUA.Nnodes sparse mapping matrix, where m=numel(x).
%
%                 Rows corresponding to locations outside of the FE mesh are identically
%                 zero. See the important note on this below.
%
%   Inside      : m x 1 logical array. Inside(i) is true if the location (x(i),y(i)) was
%                 found within the FE mesh, and false otherwise.
%
%   ID          : m x 1 array with the element number of the element containing the
%                 location (x(i),y(i)). NaN wherever Inside is false.
%
%   MUA         : the mesh structure, with the field MUA.TR populated. The underlying
%                 triangulation is only created if MUA.TR is empty or absent, so by
%                 collecting this output the (potentially repeated) cost of creating the
%                 triangulation is avoided on subsequent calls.
%
%
%% IMPORTANT: rows outside of the mesh
%
% O is returned as a m x Nnodes matrix, i.e. with one row for every location in (x,y),
% so that the rows of O are directly index-aligned with any vector of measurements.
% However, the rows corresponding to locations outside of the FE mesh are identically
% zero, and therefore
%
%   O * F
%
% returns 0 at those locations. A zero is indistinguishable from a genuine value of
% zero, and if used directly in a misfit term it would produce a spurious residual.
%
% Always use the Inside mask, for example:
%
%   Res = O(Inside,:)*F - Meas(Inside) ;
%
% or equivalently discard the outside locations up front.
%
%
%% Method
%
% 1) A triangulation over the corner nodes of the FE mesh is created using
%    CreateFEmeshTriRep (or reused from MUA.TR).
%
% 2) The enclosing element, and the barycentric coordinates within that element, are
%    found for each location using the MATLAB pointLocation method.
%
% 3) The values of the MUA.nod form functions at those barycentric coordinates are
%    evaluated using sr_shape_fun.
%
% 4) The form-function values are scattered into the rows of a sparse matrix, with the
%    column indices given by the global node numbers of the enclosing element,
%    MUA.connectivity(ID,:).
%
% Because all elements are straight-sided, the map from physical to reference
% coordinates is affine, and hence the barycentric coordinates obtained from the
% corner-node triangulation are exactly the natural coordinates of the 6- and 10-node
% elements as well. This is what makes step 3 valid for all element types.
%
% The node ordering is consistent: CreateFEmeshTriRep builds the triangulation from
% local nodes (1,2,3), (1,3,5) and (1,4,7) for nod=3, 6 and 10 respectively, and these
% are precisely the local nodes associated with the barycentric coordinates (c1,c2,c3)
% in sr_shape_fun. No permutation of the barycentric coordinates is required.
%
%
%% Note on locations lying exactly on element edges
%
% MATLAB's pointLocation returns NaN for points that lie exactly on an edge shared by
% two elements, because the point then belongs to more than one triangle. Such points
% are here re-located by nudging them by a tiny amount in each of the four diagonal
% directions in turn, until a parent element is found. The nudge is only used to
% disambiguate the element. The barycentric coordinates are subsequently evaluated at
% the original, unshifted, location using cartesianToBarycentric, and hence no
% interpolation error whatsoever is introduced by the nudge.
%
%
%% Properties of O, and some notes on use
%
% o O has exactly MUA.nod non-zero entries per inside row, i.e. nnz(O)=MUA.nod*nInside.
%
% o The rows of O sum to unity (partition of unity of the form functions). This is
%   checked internally and a warning issued if violated.
%
% o Linear fields are reproduced exactly, i.e. O*MUA.coordinates(:,1)=x for all inside
%   locations, for all element types. Quadratic fields are reproduced exactly for
%   nod>=6, and cubic fields for nod=10.
%
% o O depends on the mesh and on the (x,y) locations only, and not on any nodal values.
%   It should therefore be built once and reused, and only rebuilt following remeshing.
%
% o For a misfit term of the form
%
%       J = (O*F-d)' * inv(Sigma) * (O*F-d) / 2
%
%   the gradient with respect to the nodal values F is
%
%       dJ/dF = O' * inv(Sigma) * (O*F-d)
%
%   and the (exact) Hessian is
%
%       d2J/dF2 = O' * inv(Sigma) * O
%
%   Note that O' scatters the residuals at the measurement locations back onto the
%   nodes of the mesh. Note also that this is a pointwise misfit and not a continuous
%   L2 misfit, and consequently no mass matrix enters.
%
% o For a vector field, apply O separately to each component, or form the block matrix
%   blkdiag(O,O) acting on [u;v].
%
%
%% Example
%
%   [O,Inside]=BuildNode2DataMap(CtrlVar,MUA,Meas.x,Meas.y) ;
%   Res=O(Inside,:)*F.b-Meas.b(Inside) ;
%
% Consistency check, a linear field must be reproduced exactly:
%
%   [O,Inside]=BuildNode2DataMap(CtrlVar,MUA,xTest,yTest) ;
%   norm(O(Inside,:)*MUA.coordinates(:,1)-xTest(Inside))
%
%
%% See also
%
% CreateFEmeshTriRep, sr_shape_fun, MapNodalVariablesFromMesh1ToMesh2UsingShapeAndScattered
%
%%

narginchk(4,4)
nargoutchk(0,4)

%% Input checks and initialisation

x=x(:) ; y=y(:) ;

if numel(x)~=numel(y)
    error("BuildNode2DataMap:InconsistentInputs", ...
        "x and y must have the same number of elements (numel(x)=%i, numel(y)=%i).",numel(x),numel(y))
end

m=numel(x) ;
n=MUA.Nnodes ;

% Defaults, returned unchanged if there is nothing to do.
O=sparse(m,n) ;
Inside=false(m,1) ;
ID=nan(m,1) ;

if m==0 || MUA.Nele==0
    return
end

if isstruct(CtrlVar) && isfield(CtrlVar,"InfoLevel") && ~isempty(CtrlVar.InfoLevel)
    InfoLevel=CtrlVar.InfoLevel ;
else
    InfoLevel=0 ;
end

%% 1) The triangulation over the corner nodes of the FE mesh
%
% Only created if not already available. Note that the triangulation contains all
% MUA.Nnodes points, but only the corner nodes appear in the connectivity list.

if ~isfield(MUA,"TR") || isempty(MUA.TR)
    MUA.TR=CreateFEmeshTriRep(MUA.connectivity,MUA.coordinates) ;
end

TR=MUA.TR ;

%% 2) A length scale for the edge-disambiguation nudge
%
% The nudge must be much smaller than the smallest element, but well above eps in
% absolute terms given that coordinates are typically in metres.

if isfield(MUA,"EleAreas") && ~isempty(MUA.EleAreas)
    EleAreas=MUA.EleAreas ;
else
    EleAreas=TriAreaFE(MUA.coordinates,MUA.connectivity) ;
end

tol=1e-5*sqrt(2*min(abs(EleAreas))) ;

%% 3) Locate the enclosing elements and the barycentric coordinates

B=nan(m,3) ;

isFinite=isfinite(x) & isfinite(y) ;   % guard against NaN/Inf in the input locations

[ID(isFinite),B(isFinite,:)]=pointLocation(TR,[x(isFinite) y(isFinite)]) ;

% Deal with locations falling exactly on an element edge, for which pointLocation
% returns NaN. Nudge in each of the four diagonal directions in turn.

Missing=isFinite & isnan(ID) ;

if any(Missing)

    Shifts=tol*[ 1 1 ; -1 -1 ; -1 1 ; 1 -1 ] ;

    for iShift=1:size(Shifts,1)

        I=find(Missing) ;

        if isempty(I)
            break
        end

        IDtest=pointLocation(TR,[x(I)+Shifts(iShift,1)  y(I)+Shifts(iShift,2)]) ;

        Found=~isnan(IDtest) ;

        if any(Found)

            J=I(Found) ;
            ID(J)=IDtest(Found) ;

            % The barycentric coordinates are evaluated at the original, unshifted,
            % locations. The nudge is only used to select the parent element.
            B(J,:)=cartesianToBarycentric(TR,IDtest(Found),[x(J) y(J)]) ;

            Missing(J)=false ;

        end
    end
end

Inside=~isnan(ID) ;
nInside=sum(Inside) ;

%% 4) Form functions at those barycentric coordinates, and assembly of O

if nInside>0

    sfun=sr_shape_fun(B(Inside,:),MUA.nod) ;      % nInside x MUA.nod

    iRow=repmat(find(Inside),1,MUA.nod) ;         % nInside x MUA.nod, global row (measurement) index
    jCol=MUA.connectivity(ID(Inside),:) ;         % nInside x MUA.nod, global column (node) index

    O=sparse(iRow(:),jCol(:),sfun(:),m,n) ;

end

%% 5) Checks and information

if nInside>0

    RowSum=full(sum(O,2)) ;
    PartitionOfUnityError=max(abs(RowSum(Inside)-1)) ;

    if PartitionOfUnityError>1e-10
        warning("BuildNode2DataMap:PartitionOfUnity", ...
            "Rows of the mapping matrix do not sum to unity (max deviation=%g). \n"+ ...
            "This suggests an inconsistency between the barycentric coordinates and the form functions.", ...
            PartitionOfUnityError)
    end
end

if InfoLevel>=1 && nInside<m
    fprintf("BuildNode2DataMap: %i out of %i locations (%4.1f%%) are outside of the FE mesh and have been masked out. \n", ...
        m-nInside,m,100*(m-nInside)/m)
end

if InfoLevel>=10
    fprintf("BuildNode2DataMap: O is %i x %i with %i non-zero entries (%i per inside row). \n", ...
        m,n,nnz(O),MUA.nod)
end

end
