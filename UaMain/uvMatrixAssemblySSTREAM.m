




function [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F,BCs)


narginchk(4,4)

if nargin<5 ; Cache=[] ; end

if ~isfield(CtrlVar.uvAssembly,"SkipSymmetrisation") ; CtrlVar.uvAssembly.SkipSymmetrisation=false ; end

nargoutchk(1,4)


%%
% Assembles the Newton-Raphson system for the uv solve.
%
% Kuv : tangent matrix, where Kuv is the directional derivative of Ruv in the direction (Delta u, Delta v)
%
% The vector Ruv is basically $$<F,\phi_i>$$ 
%
% where F is the forward model.
%
% and matrix Kuv is the directional derivative of Ruv.
%
%
% Ruv=Tint-Fext;
%
% Tint   : internal nodal forces
%
% Fint   : external nodal forces
%
% Input variables:
%
% MUA  : a structure that contains data related the the FE mesh
%        For example:
%
% MUA.coordinates   : is as nNodes times 2 array with the x,y coordinates of the nodes
%
% MUA.connectivity  : contains the mesh connectivity, for a 3-nod element this would be a nEle times 3 array
%
% MUA.Nele          : Number of elements in the mesh
%
% MUA.Nnodes        : Number of nodes in the mesh
%
%
% F                 : a structure that contains all the field variables, for example:
%
% F.ub              : are the x velocity components (ub)
%
% F.vb              : the y velocity components (vb)
%
% F.h               : is the ice thickness (h)
%
%
% The Newton-Raphson system is:
%
% $$K_{uv} \, \Delta x_i = -R_{uv} $$ 
%
% $$x_{i+1}= x_{i}+ \Delta x_i $$
%
% where $i$ is the Newton-Raphson iteration number
%
% In the case of Dirichlet boundary conditions, the constraints are introduced using Lagrange variables. 
% 
% The constraints are
% 
% Luv [u,v]=Lrhsuv
% 
% where Luv is the multi-linear
% 
% (I really should rename this to something like Aeq [u;v] = b )
%
%
% I need to solve
%
% [Kxu Kxv Luv'] [du]        =   -Ru  - Luv' lambdauv
% [Kyu Kyv     ] [dv]        =   -Rv 
% [  Luv      0] [dlambdauv]     Lrhsuv-Luv [u ; v ]
%
% All matrices are Nnodes x Nnodes, apart from:
% Luv is #uv constraints x 2 Nnodes
%
%
% I write the system as
% [Kuv  Luv^T ]  [duv]  =  [ -R(uv) - Luv^T l]
% [Luv   0    ]  [dl]      [cuv-Luv uv]
%
%
%
% If the forward model is
%
% $$F(\mathbf{v}(x,y))=0 $$
%
% Then
%
% $$R_j = \langle F , \phi_i \rangle $$
%
% and
%
% $$K_{ij}=\delta_{u}R_i[\phi_j]$$
%
% and same for the $v$ derivatives.
%
% FE-formulation, note how I use
%
% $$g\, \mathcal{G} \,  (\rho h -\rho_o H^{+}) \, \partial_y B =g\, \mathcal{G} \,  (\rho h -\rho_o H^{+}) \, \partial_y b $$
%
%
% $$ R^x_i=\left \langle  h \eta ( 4 \partial_x u + 2 \partial_y v) | \partial_x \phi_i \right \rangle
%     +\langle   h \eta (\partial_y u + \partial_x v)  | \partial_y \phi_i \rangle
%    + \langle t_x | \phi_i \rangle
%    - \left \langle \frac{1}{2} g \cos(\alpha) \,  (\rho h^2 -  \rho_o d^2)  \big\vert \partial_x \phi_i \right \rangle
%    + \langle g\, \mathcal{G} \, (\rho h -\rho_o H^{+}) \partial_x B | \phi_i \rangle  - \langle \rho g \sin(\alpha) \, h  | \phi_i \rangle   =0
% $$
%
% $$ R^y_i=\langle  h \eta ( 4 \partial_y v + 2 \partial_x u) | \partial_y \phi_i \rangle
%     +\langle   h \eta (\partial_x v + \partial_y u)  | \partial_x \phi_i \rangle
%    + \langle t_y | \phi_i \rangle
%    - \left \langle \frac{1}{2} g \cos(\alpha) \, (\rho h^2 -  \rho_o d^2) | \partial_y \phi_i \right \rangle
%    + \langle g\, \mathcal{G} \, (\rho h -\rho_o H^{+}) \partial_y B | \phi_i \rangle=0
% $$
%
% $$
% \left ( \begin{array}{cc}
%    \delta_u R^x & \delta_v R^x \\
%    \delta_u R^v  & \delta_v R^y 
% \end{array} \right )
% $$
%
% For example using Weertman sliding law:
%
% $$ R^x_i=\left \langle  h \eta ( 4 \partial_x u + 2 \partial_y v) | \partial_x \phi_i \right \rangle
%     +\langle   h \eta (\partial_y u + \partial_x v)  | \partial_y \phi_i \rangle
%    + \langle \mathcal{G} \, \beta^2 \, u | \phi_i \rangle
%    - \left \langle \frac{1}{2} g \cos(\alpha) \,  (\rho h^2 -  \rho_o d^2)  \big\vert \partial_x \phi_i \right \rangle
%    + \langle g\, \mathcal{G} \, (\rho h -\rho_o H^{+}) \partial_x B | \phi_i \rangle  - \langle \rho g \sin(\alpha) \, h  | \phi_i \rangle   =0
% $$
%
% $$ R^y_i=\langle  h \eta ( 4 \partial_y v + 2 \partial_x u) | \partial_y \phi_i \rangle
%     +\langle   h \eta (\partial_x v + \partial_y u)  | \partial_x \phi_i \rangle
%    + \langle \mathcal{G} \, \beta^2 \, v | \phi_i \rangle
%    - \left \langle \frac{1}{2} g \cos(\alpha) \, (\rho h^2 -  \rho_o d^2) | \partial_y \phi_i \right \rangle
%    + \langle g\, \mathcal{G} \, (\rho h -\rho_o H^{+}) \partial_y B | \phi_i \rangle=0
% $$
%
% where
%
% $$\beta^2= (C+C_0)^{-1/m} \, (u^2+v^2+v_0^2)^{(1/m-1)/2} $$
%
%




ZeroFields=CtrlVar.uvAssembly.ZeroFields;
Ronly=CtrlVar.uvMatrixAssembly.Ronly;

if Ronly
    Kuv=[];
end

if ZeroFields
    F.ub=F.ub*0;
    F.vb=F.vb*0;
    F.ub(BCs.ubFixedNode)=BCs.ubFixedValue; F.vb(BCs.vbFixedNode)=BCs.vbFixedValue;
end

if ~CtrlVar.IncludeMelangeModelPhysics
    uoint=[];
    voint=[];
    Coint=[];
    moint=[];
    uaint=[];
    vaint=[];
    Caint=[];
    maint=[];
end


% calculates the tangent matrix (K) and right-hand side (-R) in a vectorized form

%	Ruv=Tuv-Fuv;
%   Fuv are the `external forces', i.e. right-hand side of the original system
%   Tuv are the `internal forces'. The equation is considered solved once internal and external forces are
%   equal to within a given tolerance

if any(F.h<0) ; warning('MATLAB:KRTF:hnegative',' h negative ') ; end
if any(F.C<0) ; warning('MATLAB:KRTF:Cnegative',' C negative ') ; end


if CtrlVar.TestAdjointFiniteDifferenceType=="complex step differentiation"
    CtrlVar.TestForRealValues=false;
end



if CtrlVar.TestForRealValues
    if ~isreal(F.C) ; save TestSave ; error('KRTF: C not real ') ; end
end

if any(isnan(F.ub)) ; save TestSave ; error('uvMatrixAssembly:NaN','NaN in F.ub. Variables saved in TestSave.mat') ; end
if any(isnan(F.vb)) ; save TestSave ; error('uvMatrixAssembly:NaN','NaN in F.vb. Variables saved in TestSave.mat') ; end


  

g=F.g;





%Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
ndim=2; dof=2; neq=dof*MUA.Nnodes;
neqx=MUA.Nnodes ;

%[b,s]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar);

hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
snod=reshape(F.s(MUA.connectivity,1),MUA.Nele,MUA.nod);

ubnod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vbnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);


% H=F.S-F.B;


if CtrlVar.IncludeMelangeModelPhysics

    uonod=reshape(F.uo(MUA.connectivity,1),MUA.Nele,MUA.nod);
    vonod=reshape(F.vo(MUA.connectivity,1),MUA.Nele,MUA.nod);

    uanod=reshape(F.ua(MUA.connectivity,1),MUA.Nele,MUA.nod);
    vanod=reshape(F.va(MUA.connectivity,1),MUA.Nele,MUA.nod);

end



Cnod=reshape(F.C(MUA.connectivity,1),MUA.Nele,MUA.nod);
mnod=reshape(F.m(MUA.connectivity,1),MUA.Nele,MUA.nod);

if ~isempty(F.q)
    qnod=reshape(F.q(MUA.connectivity,1),MUA.Nele,MUA.nod);
end

if ~isempty(F.V0)
    V0nod=reshape(F.V0(MUA.connectivity,1),MUA.Nele,MUA.nod);
end


if ~isempty(F.muk)
    muknod=reshape(F.muk(MUA.connectivity,1),MUA.Nele,MUA.nod);
end


if CtrlVar.IncludeMelangeModelPhysics
    Conod=reshape(F.Co(MUA.connectivity,1),MUA.Nele,MUA.nod);
    monod=reshape(F.mo(MUA.connectivity,1),MUA.Nele,MUA.nod);


    Canod=reshape(F.Ca(MUA.connectivity,1),MUA.Nele,MUA.nod);
    manod=reshape(F.ma(MUA.connectivity,1),MUA.Nele,MUA.nod);
end
%end



AGlennod=reshape(F.AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
nnod=reshape(F.n(MUA.connectivity,1),MUA.Nele,MUA.nod);




Snod=reshape(F.S(MUA.connectivity,1),MUA.Nele,MUA.nod);
Bnod=reshape(F.B(MUA.connectivity,1),MUA.Nele,MUA.nod);
rhonod=reshape(F.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);

ca=cos(F.alpha); sa=sin(F.alpha);





if ~Ronly
    d1d1=zeros(MUA.Nele,MUA.nod,MUA.nod); d2d2=zeros(MUA.Nele,MUA.nod,MUA.nod);  d1d2=zeros(MUA.Nele,MUA.nod,MUA.nod); d2d1=zeros(MUA.Nele,MUA.nod,MUA.nod);
end

Tx=zeros(MUA.Nele,MUA.nod);  Ty=zeros(MUA.Nele,MUA.nod); Fx=zeros(MUA.Nele,MUA.nod);  Fy=zeros(MUA.Nele,MUA.nod);




for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1

    Deriv=MUA.Deriv(:,:,:,Iint);  % Deriv at integration points
    detJ=MUA.DetJ(:,Iint);

    Dx=reshape(Deriv(:,1,:),MUA.Nele,MUA.nod);
    Dy=reshape(Deriv(:,2,:),MUA.Nele,MUA.nod);

    % values at this integration point

    uint=ubnod*fun;
    vint=vbnod*fun;

    if CtrlVar.IncludeMelangeModelPhysics

        uoint=uonod*fun;
        voint=vonod*fun;

        uaint=uanod*fun;
        vaint=vanod*fun;

    end

    Cint=Cnod*fun;
    Cint(Cint<CtrlVar.Cmin)=CtrlVar.Cmin; % for higher order elements it is possible that Cint is less than any of the nodal values
    mint=mnod*fun;

    if ~isempty(F.q)
        qint=qnod*fun;
    else
        qint=[];
    end

    if ~isempty(F.muk)
        mukint=muknod*fun;
    else
        mukint=[];
    end

    if ~isempty(F.V0)
        V0int=V0nod*fun;
    else
        V0int=[];
    end

    if CtrlVar.IncludeMelangeModelPhysics
        Coint=Conod*fun;
        moint=monod*fun;

        Caint=Canod*fun;
        maint=manod*fun;
    end

    AGlenint=AGlennod*fun;
    AGlenint(AGlenint<CtrlVar.AGlenmin)=CtrlVar.AGlenmin;
    nint=nnod*fun;

    sint=snod*fun;
    Bint=Bnod*fun;
    Sint=Snod*fun;
    Hint=Sint-Bint;

    if CtrlVar.Calculate.Geometry=="bh-FROM-sBS"

        bint=Bint ;     %#ok<NASGU> % ~OK, except when grounded
        hint=sint-Bint; % OK

    else    % CtrlVar.Calculate.Geometry="bs-FROM-hBS" ;

        hint=hnod*fun;
        bint=sint-hint; %#ok<NASGU>

    end

    rhoint=rhonod*fun;

    %% evaluating dint, hfint, Heint and deltaint at integration points

    hfint=F.rhow*Hint./rhoint;                                   % this is linear, so fine to evaluate at int in this manner
    Heint = HeavisideApprox(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);  % important to calculate Heint and deltaint in a consistent manner
    HEint = HeavisideApprox(CtrlVar.kH,hfint-hint,CtrlVar.Hh0);

    deltaint=DiracDelta(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);      % i.e. deltaint must be the exact derivative of Heint

    Hposint = HeavisideApprox(CtrlVar.kH,Hint,CtrlVar.Hh0).*Hint;

    dint=HEint.*rhoint.*hint/F.rhow + Heint.*Hposint ;  % definition of d, applied directly at integration points

    % derivatives at this integration point for all elements

    dsdx=sum(Dx.*snod,2);   dsdy=sum(Dy.*snod,2);
    dhdx=sum(Dx.*hnod,2);   dhdy=sum(Dy.*hnod,2);
    dBdx=sum(Dx.*Bnod,2);   dBdy=sum(Dy.*Bnod,2);

    exx=sum(Dx.*ubnod,2);
    eyy=sum(Dy.*vbnod,2);
    exy=0.5*(sum(Dx.*vbnod,2)+sum(Dy.*ubnod,2));

    [taux,tauy,dtauxdu,dtauxdv,dtauydu,dtauydv] = ...
        BasalDrag(CtrlVar,MUA,Heint,deltaint,hint,Bint,Hint,rhoint,F.rhow,uint,vint,Cint,mint,uoint,voint,Coint,moint,uaint,vaint,Caint,maint,qint,g,mukint,V0int);
    [etaint,Eint]=EffectiveViscositySSTREAM(CtrlVar,AGlenint,nint,exx,eyy,exy);

    if CtrlVar.Calculate.Geometry=="bh-FROM-sBS"
        dbdx=dBdx;  % only OK if grounded
        dbdy=dBdy;
    else
        dbdx=dsdx-dhdx;
        dbdy=dsdy-dhdy;
    end

    detJw=detJ*MUA.weights(Iint);

    %% tangent matrix

    if ~Ronly

        he=hint.*etaint.*detJw;

        tuu=dtauxdu.*detJw;
        tvv=dtauydv.*detJw;
        tuv=dtauxdv.*detJw;
        tvu=dtauydu.*detJw;

        a1=hint.*(4*exx+2*eyy);
        a2=2*hint.*exy;
        a3=hint.*(4*eyy+2*exx);

        % j-dependent, Nele x nod
        Deu=Eint.*((2*exx+eyy).*Dx+exy.*Dy);
        Dev=Eint.*((2*eyy+exx).*Dy+exy.*Dx);

        % i-dependent, Nele x nod, quadrature weight folded in
        GuW=(a1.*Dx+a2.*Dy).*detJw;
        GvW=(a3.*Dy+a2.*Dx).*detJw;

        for Inod=1:MUA.nod

            PP=fun(Inod)*fun.';        % 1 x nod

            dxi=Dx(:,Inod);   dyi=Dy(:,Inod);
            gu=GuW(:,Inod);   gv=GvW(:,Inod);

            d1d1(:,:,Inod)=d1d1(:,:,Inod) + 4*he.*dxi.*Dx + he.*dyi.*Dy + tuu.*PP + gu.*Deu;
            d2d2(:,:,Inod)=d2d2(:,:,Inod) + 4*he.*dyi.*Dy + he.*dxi.*Dx + tvv.*PP + gv.*Dev;
            d1d2(:,:,Inod)=d1d2(:,:,Inod) + he.*(2*dxi.*Dy+dyi.*Dx)     + tuv.*PP + gu.*Dev;
            d2d1(:,:,Inod)=d2d1(:,:,Inod) + he.*(2*dyi.*Dx+dxi.*Dy)     + tvu.*PP + gv.*Deu;

        end

    end

    %% residual (Tint) and external forces (Fext), vectorised over Inod

    qx =(-F.g*(rhoint.*hint-F.rhow*dint).*dbdx*ca + rhoint.*F.g.*hint.*sa).*detJw;   % Fext
    qy =(-F.g*(rhoint.*hint-F.rhow*dint).*dbdy*ca).*detJw;                           % Fext
    p2 =(0.5*F.g.*ca*(rhoint.*hint.^2-F.rhow.*dint.^2)).*detJw;                      % Fext
    p3x=(hint.*etaint.*(4*exx+2*eyy)).*detJw;                                        % Tint
    p3y=(hint.*etaint.*(4*eyy+2*exx)).*detJw;                                        % Tint
    p4 =(hint.*etaint.*2.*exy).*detJw;                                               % Tint
    tx =taux.*detJw;                                                                 % Tint
    ty =tauy.*detJw;                                                                 % Tint

    Tx=Tx + p3x.*Dx + p4.*Dy + tx.*fun.';
    Fx=Fx + qx.*fun.' + p2.*Dx;

    Ty=Ty + p3y.*Dy + p4.*Dx + ty.*fun.';
    Fy=Fy + qy.*fun.' + p2.*Dy;

end





if ~isfield(MUA,"uvAssemblyPattern") || isempty(MUA.uvAssemblyPattern)
    CtrlVar.MUA.AssemblyPattern.uv=true;
    CtrlVar.MUA.AssemblyPattern.uvh=false;
    MUA.uvAssemblyPattern=AssemblyPatternCache(CtrlVar,MUA);
end

Cache= MUA.uvAssemblyPattern;

iR=Cache.iR;
One=ones(1,1,"uint32");

Tval=zeros(MUA.nod*MUA.Nele*2,1);
Fval=zeros(MUA.nod*MUA.Nele*2,1);
istak=0;

for Inod=1:MUA.nod

    Tval(istak+1:istak+MUA.Nele)=Tx(:,Inod);
    Fval(istak+1:istak+MUA.Nele)=Fx(:,Inod);
    istak=istak+MUA.Nele;

    Tval(istak+1:istak+MUA.Nele)=Ty(:,Inod);
    Fval(istak+1:istak+MUA.Nele)=Fy(:,Inod);
    istak=istak+MUA.Nele;

end

Tint=sparse(iR,One,Tval,neq,1);
Fext=sparse(iR,One,Fval,neq,1);

Ruv=Tint-Fext;

if ~Ronly

    Xval=[d1d1(:) ; d2d2(:) ; d1d2(:) ; d2d1(:)];

    vs=accumarray(Cache.map,Xval,[Cache.nk 1]);

    if CtrlVar.uvAssembly.SkipSymmetrisation
        % diagnostic path only
    else
        vs=0.5*(vs+vs(Cache.perm));
    end

    Kuv=sparse(Cache.i0,Cache.j0,vs,neq,neq);

end




end
