
function [UserVar,M,D,rc,rd]=GeneralizedMassAndStiffnessAssembly(UserVar,CtrlVar,MUA,a,b,c,d)

%%
%
% Assembles:
%
% $$M_{pq} = \langle a\phi_p   , \phi_q \rangle  $$
%
% $$D_{pq} = \langle b \partial_x \phi_p , \partial_x \phi_q \rangle + \langle  b \partial_y \phi_p , \partial_y \phi_q \rangle  $$
%
% $$rc_{p} = \langle c , \phi_p \rangle $$
%
% $$rd_{p} = \langle \nabla d , \phi_p \rangle $$
%
% where a, b, c and d are nodal variables.
%
% One can think of M and D as generalized mass and stiffness matrices. The standard mass and stiffness matrices are
% obtained for a=1 and b=1, and typically arise when one deals with differential equations with constant coefficients.
%
%
%%



narginchk(7,7)


ndim=2; dof=1; neq=dof*MUA.Nnodes;

a=a+zeros(MUA.Nnodes,1);
b=b+zeros(MUA.Nnodes,1);
c=c+zeros(MUA.Nnodes,1);
d=d+zeros(MUA.Nnodes,1);

anod=reshape(a(MUA.connectivity,1),MUA.Nele,MUA.nod);
bnod=reshape(b(MUA.connectivity,1),MUA.Nele,MUA.nod);
cnod=reshape(c(MUA.connectivity,1),MUA.Nele,MUA.nod);
dnod=reshape(d(MUA.connectivity,1),MUA.Nele,MUA.nod);


ElementMatrixM=zeros(MUA.Nele,MUA.nod,MUA.nod);
ElementMatrixD=zeros(MUA.Nele,MUA.nod,MUA.nod);

ElementRHSc=zeros(MUA.Nele,MUA.nod);
ElementRHSd=zeros(MUA.Nele,MUA.nod);

% vector over all elements for each integration point
for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points

    Deriv=MUA.Deriv(:,:,:,Iint);
    detJ=MUA.DetJ(:,Iint);



    % Deriv : Nele x dof x nod
    %  detJ : Nele

    % values at integration point

    aint=anod*fun;
    bint=bnod*fun;
    cint=cnod*fun;

    dddx=zeros(MUA.Nele,1);
    dddy=zeros(MUA.Nele,1);

    for Inod=1:MUA.nod

        dddx=dddx+Deriv(:,1,Inod).*dnod(:,Inod);
        dddy=dddy+Deriv(:,2,Inod).*dnod(:,Inod);

    end


    detJw=detJ*MUA.weights(Iint);

    % dt theta ( d(u1 h1)/dx    + d(v1 h1)/dy) + h1=
    %  h0+dt { (1-theta) a0+theta a1-(1-theta) (d(u0 h0)/dx+d(v0 h0)/dy}



    for Inod=1:MUA.nod


        for Jnod=1:MUA.nod

            aphiphi=aint.*fun(Jnod).*fun(Inod).*detJw;

            bdphidxdphidx=bint.*Deriv(:,1,Jnod).*Deriv(:,1,Inod).*detJw;
            bdphidydphidy=bint.*Deriv(:,2,Jnod).*Deriv(:,2,Inod).*detJw;

            ElementMatrixM(:,Inod,Jnod)=ElementMatrixM(:,Inod,Jnod)+aphiphi ;
            ElementMatrixD(:,Inod,Jnod)=ElementMatrixD(:,Inod,Jnod)+bdphidxdphidx+bdphidydphidy ;

        end

        ElementRHSc(:,Inod)=ElementRHSc(:,Inod)+...
            cint.*fun(Inod).*detJw ;


        ElementRHSd(:,Inod)=ElementRHSd(:,Inod)+...
            (dddx.*Deriv(:,1,Inod)+ dddy.*Deriv(:,2,Inod)).*detJw ;

    end
end




Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1); Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
XvalM=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
XvalD=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
istak=0;

for Inod=1:MUA.nod
    for Jnod=1:MUA.nod
        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);
        XvalM(istak+1:istak+MUA.Nele)=ElementMatrixM(:,Inod,Jnod);
        XvalD(istak+1:istak+MUA.Nele)=ElementMatrixD(:,Inod,Jnod);
        istak=istak+MUA.Nele;
    end
end

M=sparseUA(Iind,Jind,XvalM,neq,neq);

if nargout>2
    D=sparseUA(Iind,Jind,XvalD,neq,neq);
else
    D=[];
end

if nargout>3

    % assemble right-hand sides

    rc=sparseUA(neq,1);
    for Inod=1:MUA.nod
        rc=rc+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),ElementRHSc(:,Inod),neq,1);
    end


    rd=sparseUA(neq,1);
    for Inod=1:MUA.nod
        rd=rd+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),ElementRHSd(:,Inod),neq,1);
    end
else
    rc=[] ; rd=[];
end



end






















