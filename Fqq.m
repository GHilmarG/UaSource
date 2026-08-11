


function KFqq=Fqq(CtrlVar,MUA,F,uAdjoint,vAdjoint)

narginchk(5,5)

%

narginchk(5,5)

ndim=2;
nNodes=MUA.Nnodes ;


hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);
uAdjointnod=reshape(uAdjoint(MUA.connectivity,1),MUA.Nele,MUA.nod);
vAdjointnod=reshape(vAdjoint(MUA.connectivity,1),MUA.Nele,MUA.nod);
AGlennod=reshape(F.AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
nnod=reshape(F.n(MUA.connectivity,1),MUA.Nele,MUA.nod);


H=zeros(MUA.Nele,MUA.nod,MUA.nod);

for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    detJ=MUA.DetJ(:,Iint);
    Deriv=MUA.Deriv(:,:,:,Iint);

    h=hnod*fun;
    n=nnod*fun;
    A=AGlennod*fun;
    A(A<CtrlVar.AGlenmin)=CtrlVar.AGlenmin;


    exx=zeros(MUA.Nele,1); exy=zeros(MUA.Nele,1);  eyy=zeros(MUA.Nele,1);
    dlxdx=zeros(MUA.Nele,1); dlydx=zeros(MUA.Nele,1); dlxdy=zeros(MUA.Nele,1); dlydy=zeros(MUA.Nele,1);

    for Inod=1:MUA.nod

        exx=exx+Deriv(:,1,Inod).*unod(:,Inod);
        eyy=eyy+Deriv(:,2,Inod).*vnod(:,Inod);
        exy=exy+0.5*(Deriv(:,1,Inod).*vnod(:,Inod) + Deriv(:,2,Inod).*unod(:,Inod));


        dlxdx=dlxdx+Deriv(:,1,Inod).*uAdjointnod(:,Inod);
        dlydx=dlydx+Deriv(:,1,Inod).*vAdjointnod(:,Inod);
        dlxdy=dlxdy+Deriv(:,2,Inod).*uAdjointnod(:,Inod);
        dlydy=dlydy+Deriv(:,2,Inod).*vAdjointnod(:,Inod);

    end


    detJw=detJ*MUA.weights(Iint);

    e=real(sqrt(CtrlVar.EpsZero^2+exx.^2+eyy.^2+exx.*eyy+exy.^2));
    detadA=   real(  A.^(-1./n-1) .* e.^(1./n-1) .* 0.5* (-1./n)            ); 
    d2etadAdA=real(  A.^(-1./n-2) .* e.^(1./n-1) .* 0.5.*(1./(n.^2)+(1./n)) );
 %   Temp=d2etadAdA.*h.*((4*exx+2*eyy).*dlxdx+2*exy.*dlxdy+(4*eyy+2*exx).*dlydy+2*exy.*dlydx) .* (log(10).*A).^2; 

    Temp=-h.*( (4*exx+2*eyy).*dlxdx + 2*exy.*dlxdy + (4*eyy+2*exx).*dlydy + 2*exy.*dlydx);

    dfdx=detadA.*Temp;
    d2fdxdx=d2etadAdA.*Temp; 

    d2fdyy=log(10)^2 .* ( A.^2 .* d2fdxdx + A.*dfdx);
    for Inod=1:MUA.nod
        for Jnod=1:MUA.nod

            H(:,Inod,Jnod)=H(:,Inod,Jnod) + d2fdyy .*fun(Inod) .*fun(Jnod).*detJw;


            
            Fuu=h.*d2uu_eta_ik.*(4.*dudx+2.*dvdy).*dlxdx

        end
    end
end


Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1,'uint32');
Xval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);

istak=0;

for Inod=1:MUA.nod
 

    for Jnod=1:MUA.nod

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);
        Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);
        Xval(istak+1:istak+MUA.Nele)=H(:,Inod,Jnod);
        %Xval(istak+1:istak+MUA.Nele)=H(:,Inod,Jnod).* (log(10).*F.AGlen(MUA.connectivity(:,Inod))) .*(log(10).*F.AGlen(MUA.connectivity(:,Jnod))) ; % At nodes

        istak=istak+MUA.Nele;

    end
end


HPsiddFdAdA=sparseUA(Iind,Jind,Xval,nNodes,nNodes);

HPsiddFdAdA=(HPsiddFdAdA'+HPsiddFdAdA)/2 ;




end







