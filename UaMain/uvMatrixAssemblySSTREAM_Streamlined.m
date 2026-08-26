




function [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM_Streamlined(CtrlVar,MUA,F,BCs)

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


narginchk(4,4)
nargoutchk(1,4)



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


    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points


    Deriv=MUA.Deriv(:,:,:,Iint);  % Deriv at integration points
    detJ=MUA.DetJ(:,Iint);



    %        fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    %       [Deriv,detJ]=derivVector(coordinates,connectivity,nip,Iint);

    % Deriv : Nele x dof x nod
    %  detJ : Nele

    % values at integration this point


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
    %   end


    AGlenint=AGlennod*fun;
    AGlenint(AGlenint<CtrlVar.AGlenmin)=CtrlVar.AGlenmin;
    nint=nnod*fun;
  

    % hint=hnod*fun;
    % sint=snod*fun;
    % Bint=Bnod*fun;
    % Sint=Snod*fun;
    % bint=sint-hint;
    % Hint=Sint-Bint;


    sint=snod*fun;
    Bint=Bnod*fun;
    Sint=Snod*fun;
    Hint=Sint-Bint;

    if CtrlVar.Calculate.Geometry=="bh-FROM-sBS"

        bint=Bint ;     % ~OK, except when grounded
        hint=sint-bint; % OK

    else    % CtrlVar.Calculate.Geometry="bs-FROM-hBS" ;

        hint=hnod*fun;  %  I could put calculating bs from hBS in here
        bint=sint-hint; %

    end

    rhoint=rhonod*fun;



    %% evaluating dint, hfint, Heint and deltaint at integration points#



    % $$ d=\mathcal{H}(h_f-h) \, \rho h /\rho_w + \mathcal{H}(h-h_f) \,  H^{+} $$
    % 2024/12/28: This is the new post 2025 default. This is consistent with same terms in the uvh assembly

    hfint=F.rhow*Hint./rhoint;                                   % this is linear, so fine to evaluate at int in this manner
    Heint = HeavisideApprox(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);  % important to calculate Heint and deltaint in a consistent manner
    HEint = HeavisideApprox(CtrlVar.kH,hfint-hint,CtrlVar.Hh0);

    deltaint=DiracDelta(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);       % i.e. deltaint must be the exact derivative of Heint
  

    Hposint = HeavisideApprox(CtrlVar.kH,Hint,CtrlVar.Hh0).*Hint;

    % Here we apply the definition of d directly at integration points
    dint=HEint.*rhoint.*hint/F.rhow + Heint.*Hposint ;  % definition of d, applied directly at integration points


    % derivatives at this integration point for all elements
    dsdx=zeros(MUA.Nele,1); dhdx=zeros(MUA.Nele,1);
    dsdy=zeros(MUA.Nele,1); dhdy=zeros(MUA.Nele,1);
    dBdx=zeros(MUA.Nele,1); dBdy=zeros(MUA.Nele,1);


    exx=zeros(MUA.Nele,1);
    eyy=zeros(MUA.Nele,1);
    exy=zeros(MUA.Nele,1);


    for Inod=1:MUA.nod

        dsdx=dsdx+Deriv(:,1,Inod).*snod(:,Inod);
        dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
        dsdy=dsdy+Deriv(:,2,Inod).*snod(:,Inod);
        dhdy=dhdy+Deriv(:,2,Inod).*hnod(:,Inod);

        dBdx=dBdx+Deriv(:,1,Inod).*Bnod(:,Inod);
        dBdy=dBdy+Deriv(:,2,Inod).*Bnod(:,Inod);

        exx=exx+Deriv(:,1,Inod).*ubnod(:,Inod);
        eyy=eyy+Deriv(:,2,Inod).*vbnod(:,Inod);
        exy=exy+0.5*(Deriv(:,1,Inod).*vbnod(:,Inod) + Deriv(:,2,Inod).*ubnod(:,Inod));


    end



    [taux,tauy,dtauxdu,dtauxdv,dtauydu,dtauydv] = ...
        BasalDrag(CtrlVar,MUA,Heint,deltaint,hint,Bint,Hint,rhoint,F.rhow,uint,vint,Cint,mint,uoint,voint,Coint,moint,uaint,vaint,Caint,maint,qint,g,mukint,V0int);
    [etaint,Eint]=EffectiveViscositySSTREAM(CtrlVar,AGlenint,nint,exx,eyy,exy);

    %%
    % Do we calculate : b and h from sBS, or b and s from hBS ?
    %
    %
    % The typical situation is to consider hBS as independent variables, and to calculate s (upper glacier surface) and b (lower
    % glacier surface) from h (ice thickness), B (bedrock) and S (ocean surface). And this is the option in all diagnostic uv solves.
    % In this case 
    % CtrlVar.Calculate.Geometry="bs-FROM-hBS" ;
    %
    % However, when doing an inversion with respect to B where s is based on measurements, the approach is to consider s, B and S
    % as independent and b and h as functions of s, B and S (and densities). In this case
    % CtrlVar.MapOldToNew.Transient.Geometry="bh-FROM-sBS" ; 
    %
    
  
    if CtrlVar.Calculate.Geometry=="bh-FROM-sBS"
        dbdx=dBdx;  % only OK if grounded
        dbdy=dBdy;
    else
        dbdx=dsdx-dhdx;   %  I could put calculating bs from hBS in here
        dbdy=dsdy-dhdy;
    end
    %%
  

    detJw=detJ*MUA.weights(Iint);


    for Inod=1:MUA.nod
        if ~Ronly
            for Jnod=1:MUA.nod


                d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)...
                    +(4*hint.*etaint.*Deriv(:,1,Inod).*Deriv(:,1,Jnod)...
                    +hint.*etaint.*Deriv(:,2,Inod).*Deriv(:,2,Jnod)...
                    +dtauxdu.*fun(Jnod).*fun(Inod)...
                    ).*detJw;


                d2d2(:,Inod,Jnod)=d2d2(:,Inod,Jnod)...
                    +(4*hint.*etaint.*Deriv(:,2,Inod).*Deriv(:,2,Jnod)...
                    +hint.*etaint.*Deriv(:,1,Inod).*Deriv(:,1,Jnod)...
                    +dtauydv.*fun(Jnod).*fun(Inod)...
                    ).*detJw ;



                d1d2(:,Inod,Jnod)=d1d2(:,Inod,Jnod)...
                    +(etaint.*hint.*(2*Deriv(:,1,Inod).*Deriv(:,2,Jnod)+Deriv(:,2,Inod).*Deriv(:,1,Jnod))...
                    + +dtauxdv.*fun(Jnod).*fun(Inod)...
                    ).*detJw;


                d2d1(:,Inod,Jnod)=d2d1(:,Inod,Jnod)...
                    +(etaint.*hint.*(2*Deriv(:,2,Inod).*Deriv(:,1,Jnod)+Deriv(:,1,Inod).*Deriv(:,2,Jnod))...
                    +dtauydu.*fun(Jnod).*fun(Inod)...
                    ).*detJw;

                %                dxu=E (2 exx+eyy)
                %                dyu=E exy
                %                dyv=E (2 eyy + exx )
                %                dxv=E exy = dyu

                Deu=Eint.*((2*exx+eyy).*Deriv(:,1,Jnod)+exy.*Deriv(:,2,Jnod));
                Dev=Eint.*((2*eyy+exx).*Deriv(:,2,Jnod)+exy.*Deriv(:,1,Jnod));

                % E11=h Deu (4 p_x u + 2 p_y v)   + h Deu  ( p_x v + p_y u) p_y N_p

                E11=  hint.*(4.*exx+2.*eyy).*Deu.*Deriv(:,1,Inod)...
                    +2*hint.*exy.*Deu.*Deriv(:,2,Inod);


                E12=  hint.*(4.*exx+2.*eyy).*Dev.*Deriv(:,1,Inod)...
                    +2*hint.*exy.*Dev.*Deriv(:,2,Inod);



                E22=  hint.*(4.*eyy+2.*exx).*Dev.*Deriv(:,2,Inod)...
                    +2*hint.*exy.*Dev.*Deriv(:,1,Inod);


                E21= hint.*(4.*eyy+2.*exx).*Deu.*Deriv(:,2,Inod)...
                    +2*hint.*exy.*Deu.*Deriv(:,1,Inod);



                d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+E11.*detJw;
                d2d2(:,Inod,Jnod)=d2d2(:,Inod,Jnod)+E22.*detJw;
                d1d2(:,Inod,Jnod)=d1d2(:,Inod,Jnod)+E12.*detJw;
                d2d1(:,Inod,Jnod)=d2d1(:,Inod,Jnod)+E21.*detJw;

            end
        end


        % R=Tint-Fext 
        t1=-F.g*    (rhoint.*hint-F.rhow*dint).*dbdx.*fun(Inod)*ca+ rhoint.*F.g.*hint.*sa.*fun(Inod);  % Fext
        t2=0.5*F.g.*ca*(rhoint.*hint.^2-F.rhow.*dint.^2).*Deriv(:,1,Inod);                             % Fext
        t3=hint.*etaint.*(4*exx+2*eyy).*Deriv(:,1,Inod);                                               % Tint
        t4=hint.*etaint.*2.*exy.*Deriv(:,2,Inod);                                                      % Tint
        t5=taux.*fun(Inod);                                                                            % Tint

        Tx(:,Inod)=Tx(:,Inod)+(t3+t4+t5).*detJw;
        Fx(:,Inod)=Fx(:,Inod)+(t1+t2).*detJw;



        t1=-F.g*(rhoint.*hint-F.rhow*dint).*dbdy.*fun(Inod)*ca;
        t2=0.5*ca*F.g.*(rhoint.*hint.^2-F.rhow.*dint.^2).*Deriv(:,2,Inod);
        t3=hint.*etaint.*(4*eyy+2*exx).*Deriv(:,2,Inod);
        t4=hint.*etaint.*2.*exy.*Deriv(:,1,Inod);
        t5=tauy.*fun(Inod);

        Ty(:,Inod)=Ty(:,Inod)+(t3+t4+t5).*detJw;
        Fy(:,Inod)=Fy(:,Inod)+(t1+t2).*detJw;





    end
end


iR=zeros(MUA.nod*MUA.Nele*2,1,"uint32");
One=ones(1,1,"uint32");
Tval=zeros(MUA.nod*MUA.Nele*2,1);
Fval=zeros(MUA.nod*MUA.Nele*2,1);
istak=0;

for Inod=1:MUA.nod


    iR(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);
    Tval(istak+1:istak+MUA.Nele)=Tx(:,Inod);
    Fval(istak+1:istak+MUA.Nele)=Fx(:,Inod);

    istak=istak+MUA.Nele;
    iR(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+neqx;
    Tval(istak+1:istak+MUA.Nele)=Ty(:,Inod);
    Fval(istak+1:istak+MUA.Nele)=Fy(:,Inod);

    istak=istak+MUA.Nele;

end
Tint=sparseUA(iR,One,Tval,neq,1);
Fext=sparseUA(iR,One,Fval,neq,1);



Ruv=Tint-Fext;

if ~Ronly


    % collect all entries in vectors, and call sparse just once


    Iind=zeros(MUA.nod*MUA.nod*MUA.Nele*4,1,'uint32'); Jind=zeros(MUA.nod*MUA.nod*MUA.Nele*4,1,'uint32');

    Xval=zeros(MUA.nod*MUA.nod*MUA.Nele*4,1);
    istak=0;

    for Inod=1:MUA.nod

        for Jnod=1:MUA.nod

            Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=d1d1(:,Inod,Jnod);
            istak=istak+MUA.Nele;

            Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+neqx; Xval(istak+1:istak+MUA.Nele)=d2d2(:,Inod,Jnod);
            istak=istak+MUA.Nele;

            Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod)+neqx; Xval(istak+1:istak+MUA.Nele)=d1d2(:,Inod,Jnod);
            istak=istak+MUA.Nele;

            Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+neqx; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=d2d1(:,Inod,Jnod);
            %Iind(istak+1:istak+Nele)=connectivity(:,Inod)+neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod); Xval(istak+1:istak+Nele)=d1d2(:,Jnod,Inod);
            istak=istak+MUA.Nele;

        end

    end


    Kuv=sparse(Iind,Jind,Xval,neq,neq);

    if CtrlVar.TestForRealValues
        Kuv=(Kuv+Kuv.')/2 ;
    else
        Kuv=(Kuv+Kuv')/2 ;
    end

end




end



