function [N,xp,info]=PairwiseConstraintsNullSpace(Aeq,beq,options)

%%
%   [N,xp,info]=PairwiseConstraintsNullSpace(Aeq,beq)
%
% Parametrises the solution set of the linear constraints
%
%   Aeq x = beq
%
% as
%
%   x = xp + N y
%
% where y is a vector of free variables, xp is a particular solution, and the columns of N are an orthonormal basis for the
% null space of Aeq.
%
% This works if each row of Aeq has at most two non-zero entries, i.e. if every constraint is either
%
%   a x_i = b               (fixed value)
%
% or
%
%   a_1 x_i + a_2 x_j = b   (pairwise tie, or e.g. a normal-velocity condition n_x u + n_y v = b).
%
% The coefficients a, a_1 and a_2 and the right-hand side b can be arbitrary (apart from being non-zero in the case of the
% coefficients), so scaled rows and ties with offsets are fine.
%
% Method: Each constraint with two entries is an affine relation x_j = s x_i + t between two dofs. These relations define a
% graph. In each connected component all dofs are affine functions of one variable, z, say of the value at the root of the
% component, x_k = s_k z + t_k. If a component contains a fixed dof, or a loop that fixes z, then z is known and the whole
% component is fixed. Otherwise the component contributes one column to N, with entries s_k/norm(s).
%
% Redundant but consistent constraints are fine. Inconsistent constraints give an error.
%
% Outputs:
%
%   N     (nx x ny) sparse matrix with orthonormal columns, N'*N=I.
%   xp    (nx x 1) particular solution, Aeq*xp=beq.
%   info  structure with some information about the constraints, including info.Residual which is the largest of
%         |Aeq*xp-beq| and |Aeq*N|.
%
% Example (u_1 fixed, u_2=u_3, and 2 u_4 - u_5 = 1):
%
%   Aeq=sparse([1 0 0 0 0 ; 0 1 -1 0 0 ; 0 0 0 2 -1]) ; beq=[3;0;1] ;
%   [N,xp]=PairwiseConstraintsNullSpace(Aeq,beq) ;
%
% See also: uvh_GaussNewton_fsolve
%%

arguments
    Aeq {mustBeNumeric}
    beq {mustBeNumeric}
    options.Tol (1,1) double = 1e-10
end

Tol=options.Tol;

Aeq=sparse(Aeq);
beq=full(double(beq(:)));
[m,nx]=size(Aeq);

if numel(beq)~=m
    error("PairwiseConstraintsNullSpace:SizeMismatch","Number of rows in Aeq (%i) and number of elements in beq (%i) must be the same.",m,numel(beq))
end

info.nFixedRows=0;
info.nTieRows=0;
info.nFixedDofs=0;
info.nFreeGroups=nx;
info.Residual=0;

if m==0
    N=speye(nx);
    xp=zeros(nx,1);
    return
end

%% Sort the constraints according to the number of non-zero entries in each row

[dof,row,val]=find(Aeq.');                 % entries are sorted by row of Aeq, and then by dof
nnzRow=accumarray(row,1,[m 1]);

if any(nnzRow>2)
    iRow=find(nnzRow>2,1);
    error("PairwiseConstraintsNullSpace:TooManyEntries","Constraint row %i has %i non-zero entries. Only rows with at most two entries are supported.",iRow,nnzRow(iRow))
end

emptyRows=find(nnzRow==0);
if any(abs(beq(emptyRows))>Tol)
    error("PairwiseConstraintsNullSpace:Inconsistent","Constraint row %i has no non-zero entries, but a non-zero right-hand side.",emptyRows(find(abs(beq(emptyRows))>Tol,1)))
end

% fixed values:  a x_j = b
isF=nnzRow(row)==1;
dofF=dof(isF);
vF=beq(row(isF))./val(isF);

% ties:  a1 x_j1 + a2 x_j2 = b
isT=nnzRow(row)==2;
dT=reshape(dof(isT),2,[]).';
aT=reshape(val(isT),2,[]).';
rT=row(isT) ; rT=rT(1:2:end);
bT=beq(rT);
j1=dT(:,1) ; j2=dT(:,2) ;
a1=aT(:,1) ; a2=aT(:,2) ;
nTie=numel(j1);

info.nFixedRows=numel(dofF);
info.nTieRows=nTie;

%% Affine relations between dofs.  x_dst = S x_src + T
%
%   forward :  x_j2 = (-a1/a2) x_j1 + b/a2
%   backward:  x_j1 = (-a2/a1) x_j2 + b/a1
%

comp=(1:nx).';          % comp(k) is the index of the root dof of the component that dof k belongs to
s=ones(nx,1);           % x_k = s_k z + t_k, where z is the value of the root dof
t=zeros(nx,1);
zc=zeros(0,1);          % extra information on z from loops that are not redundant
zv=zeros(0,1);

if nTie>0

    src=[j1;j2];
    dst=[j2;j1];
    S=[-a1./a2 ; -a2./a1];
    T=[bT./a2 ; bT./a1];

    [src,ord]=sort(src);
    dst=dst(ord) ; S=S(ord) ; T=T(ord);
    nEdges=accumarray(src,1,[nx 1]);
    ptr=[0;cumsum(nEdges)];                 % edges of dof i are ptr(i)+1:ptr(i+1)

    visited=false(nx,1);
    queue=zeros(nx,1);
    roots=find(nEdges>0);

    for r=roots.'

        if visited(r) ; continue ; end

        visited(r)=true ; queue(1)=r ; head=1 ; tail=1 ;

        while head<=tail

            i=queue(head) ; head=head+1 ;

            for e=ptr(i)+1:ptr(i+1)

                jn=dst(e);
                sn=S(e)*s(i);
                tn=S(e)*t(i)+T(e);

                if ~visited(jn)

                    visited(jn)=true ;
                    comp(jn)=r ; s(jn)=sn ; t(jn)=tn ;
                    tail=tail+1 ; queue(tail)=jn ;

                else

                    % Edge closes a loop. The two expressions for x_jn must agree:   s(jn) z + t(jn) = sn z + tn
                    ds=s(jn)-sn ;
                    dt=tn-t(jn) ;

                    if abs(ds)<=Tol*max(1,abs(sn))
                        if abs(dt)>Tol*(1+abs(tn)+abs(t(jn)))
                            error("PairwiseConstraintsNullSpace:Inconsistent","Inconsistent constraints in a loop of ties involving dof %i and dof %i.",i,jn)
                        end
                    else
                        % the loop fixes the value of z
                        zc(end+1,1)=r ;  %#ok<AGROW>
                        zv(end+1,1)=dt/ds ; %#ok<AGROW>
                    end
                end
            end
        end
    end
end

%% Which components are fixed?  Fixed values and loops both give values for z

c=[comp(dofF);zc];
v=[(vF-t(dofF))./s(dofF);zv];

zmin=nan(nx,1);
if ~isempty(c)
    zmin=accumarray(c,v,[nx 1],@min,NaN);
    zmax=accumarray(c,v,[nx 1],@max,NaN);
    iBad=find(abs(zmax-zmin)>Tol*(1+abs(zmin)),1);
    if ~isempty(iBad)
        error("PairwiseConstraintsNullSpace:Inconsistent","Inconsistent constraints: more than one value is prescribed for the group of dofs connected to dof %i.",iBad)
    end
end

isFixedComp=~isnan(zmin);          % indexed by root dof
nodeFixed=isFixedComp(comp);

%% Particular solution

xp=t;
xp(nodeFixed)=s(nodeFixed).*zmin(comp(nodeFixed))+t(nodeFixed);

%% Null-space basis. One column for each free component

isRoot=(comp==(1:nx).');
freeRoot=isRoot & ~isFixedComp;
ny=nnz(freeRoot);

yIdx=zeros(nx,1);
yIdx(freeRoot)=1:ny;

nrm=sqrt(accumarray(comp,s.^2,[nx 1]));

k=find(~nodeFixed);
N=sparse(k,yIdx(comp(k)),s(k)./nrm(comp(k)),nx,ny);

%% Check

res1=norm(Aeq*xp-beq,inf);
res2=norm(Aeq*N,inf);
info.Residual=max(res1,res2);
info.nFixedDofs=nnz(nodeFixed);
info.nFreeGroups=ny;

if info.Residual>sqrt(Tol)*(1+norm(beq,inf))
    error("PairwiseConstraintsNullSpace:CheckFailed","The constructed parametrisation does not satisfy the constraints: max(|Aeq xp-beq|)=%g and max(|Aeq N|)=%g.",res1,res2)
end

end
