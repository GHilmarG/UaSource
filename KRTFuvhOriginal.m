function [R,K,T,F]=KRTFuvhOriginal(u,v,h,S,B,u0,v0,h0,as,ab,dt,AGlen,n,C,m,coordinates,connectivity,nip,alpha,rho,rhow,g,MeshProp,CtrlVar)


if nargout==1 ; Ronly=1; else Ronly=0;end


Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
ndim=2;  neq=3*Nnodes;
neqx=Nnodes ;

theta=CtrlVar.theta;


H=S-B;
hnod=reshape(h(connectivity,1),Nele,nod);   % Nele x nod
unod=reshape(u(connectivity,1),Nele,nod);
vnod=reshape(v(connectivity,1),Nele,nod);
Cnod=reshape(C(connectivity,1),Nele,nod);
Hnod=reshape(H(connectivity,1),Nele,nod);
Snod=reshape(S(connectivity,1),Nele,nod);
Bnod=reshape(B(connectivity,1),Nele,nod);
h0nod=reshape(h0(connectivity,1),Nele,nod);
u0nod=reshape(u0(connectivity,1),Nele,nod);
v0nod=reshape(v0(connectivity,1),Nele,nod);
asnod=reshape(as(connectivity,1),Nele,nod);
abnod=reshape(ab(connectivity,1),Nele,nod);

rhog=rho*g; ca=cos(alpha); sa=sin(alpha); varrhog=g*rho*(1-rho/rhow);


[points,weights]=sample('triangle',nip,ndim);

if ~Ronly
    Kxu=zeros(Nele,nod,nod); Kxv=zeros(Nele,nod,nod);  Kxh=zeros(Nele,nod,nod);
    Kyu=zeros(Nele,nod,nod); Kyv=zeros(Nele,nod,nod);  Kyh=zeros(Nele,nod,nod);
    Khu=zeros(Nele,nod,nod); Khv=zeros(Nele,nod,nod);  Khh=zeros(Nele,nod,nod);
end

Tx=zeros(Nele,nod);  Ty=zeros(Nele,nod); Fx=zeros(Nele,nod);  Fy=zeros(Nele,nod); Th=zeros(Nele,nod);  Fh=zeros(Nele,nod);

[etaInt,~,~,exx,eyy,exy,Eint]=calcStrainRatesEtaInt(u,v,coordinates,connectivity,nip,AGlen,n,CtrlVar);

for Iint=1:nip
    
    
    fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    %[Deriv,detJ]=derivVector(coordinates,connectivity,nip,Iint);
    Deriv=MeshProp.Deriv{Iint} ; detJ=MeshProp.detJ{Iint};
    % Deriv : Nele x dof x nod
    %  detJ : Nele
    
    % values at integration this point
    hint=hnod*fun;   uint=unod*fun; vint=vnod*fun; Cint=Cnod*fun;
    h0int=h0nod*fun; u0int=u0nod*fun; v0int=v0nod*fun;
    asint=asnod*fun; abint=abnod*fun;
    Bint=Bnod*fun; Sint=Snod*fun;
    
    % calculate He and DiracDelta from integration point values of thickness
    
    hfint=rhow*(Sint-Bint)/rho;
    %Heint=hint > hfint;  % not a good idea
    kH=CtrlVar.kH;
    Heint = HeavisideApprox(kH,hint-hfint);  % very important to calculate Heint and deltaint in consistent manner
    deltaint=DiracDelta(kH,hint-hfint);      % i.e. deltaint must be the exact derivative of Heint
    
    
    
    [beta2int,Dbeta2Duuint,Dbeta2Dvvint,Dbeta2Duvint] = calcBeta2in2Dint(uint,vint,Cint,m,Heint,CtrlVar);
    etaint=etaInt(:,Iint) ;  % I could consider calculating this here
    
    
    dhdx=zeros(Nele,1); dhdy=zeros(Nele,1); dHdx=zeros(Nele,1); dHdy=zeros(Nele,1);
    dh0dx=zeros(Nele,1); dh0dy=zeros(Nele,1);
    du0dx=zeros(Nele,1); du0dy=zeros(Nele,1);
    dv0dx=zeros(Nele,1); dv0dy=zeros(Nele,1);
    
    for Inod=1:nod
        
        dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
        dhdy=dhdy+Deriv(:,2,Inod).*hnod(:,Inod);
        
        dHdx=dHdx+Deriv(:,1,Inod).*Hnod(:,Inod);
        dHdy=dHdy+Deriv(:,2,Inod).*Hnod(:,Inod);
        
        dh0dx=dh0dx+Deriv(:,1,Inod).*h0nod(:,Inod);
        dh0dy=dh0dy+Deriv(:,2,Inod).*h0nod(:,Inod);
        
        du0dx=du0dx+Deriv(:,1,Inod).*u0nod(:,Inod);
        du0dy=du0dy+Deriv(:,2,Inod).*u0nod(:,Inod);
        
        dv0dx=dv0dx+Deriv(:,1,Inod).*v0nod(:,Inod);
        dv0dy=dv0dy+Deriv(:,2,Inod).*v0nod(:,Inod);
    end
    
    
    detJw=detJ*weights(Iint);
    
    
    for Inod=1:nod
        if ~Ronly
            for Jnod=1:nod
                
                Deu=Eint(:,Iint).*((2*exx(:,Iint)+eyy(:,Iint)).*Deriv(:,1,Jnod)+exy(:,Iint).*Deriv(:,2,Jnod));
                Dev=Eint(:,Iint).*((2*eyy(:,Iint)+exx(:,Iint)).*Deriv(:,2,Jnod)+exy(:,Iint).*Deriv(:,1,Jnod));
                % E11=h Deu (4 p_x u + 2 p_y v)   + h Deu  ( p_x v + p_y u) p_y N_p
                E11=  hint.*(4.*exx(:,Iint)+2.*eyy(:,Iint)).*Deu.*Deriv(:,1,Inod)+2*hint.*exy(:,Iint).*Deu.*Deriv(:,2,Inod);
                E12=  hint.*(4.*exx(:,Iint)+2.*eyy(:,Iint)).*Dev.*Deriv(:,1,Inod)+2*hint.*exy(:,Iint).*Dev.*Deriv(:,2,Inod);
                E22=  hint.*(4.*eyy(:,Iint)+2.*exx(:,Iint)).*Dev.*Deriv(:,2,Inod)+2*hint.*exy(:,Iint).*Dev.*Deriv(:,1,Inod);
                E21=  hint.*(4.*eyy(:,Iint)+2.*exx(:,Iint)).*Deu.*Deriv(:,2,Inod)+2*hint.*exy(:,Iint).*Deu.*Deriv(:,1,Inod);
                
                
                Kxu(:,Inod,Jnod)=Kxu(:,Inod,Jnod)...
                    +(4*hint.*etaint.*Deriv(:,1,Inod).*Deriv(:,1,Jnod)...
                    +hint.*etaint.*Deriv(:,2,Inod).*Deriv(:,2,Jnod)...
                    +beta2int.*fun(Jnod).*fun(Inod)...
                    +E11...
                    +Dbeta2Duuint.*fun(Jnod).*fun(Inod)).*detJw;
                
                Kyv(:,Inod,Jnod)=Kyv(:,Inod,Jnod)...
                    +(4*hint.*etaint.*Deriv(:,2,Inod).*Deriv(:,2,Jnod)...
                    +hint.*etaint.*Deriv(:,1,Inod).*Deriv(:,1,Jnod)...
                    +beta2int.*fun(Jnod).*fun(Inod)...
                    +E22...
                    +Dbeta2Dvvint.*fun(Jnod).*fun(Inod)).*detJw ;
                
                
                
                Kxv(:,Inod,Jnod)=Kxv(:,Inod,Jnod)...
                    +(etaint.*hint.*(2*Deriv(:,1,Inod).*Deriv(:,2,Jnod)+Deriv(:,2,Inod).*Deriv(:,1,Jnod))...
                    +E12...
                    +Dbeta2Duvint.*fun(Jnod).*fun(Inod)).*detJw;    % beta derivative
                
                
                Kyu(:,Inod,Jnod)=Kyu(:,Inod,Jnod)...
                    +(etaint.*hint.*(2*Deriv(:,2,Inod).*Deriv(:,1,Jnod)+Deriv(:,1,Inod).*Deriv(:,2,Jnod))...
                    +E21...
                    +Dbeta2Duvint*fun(Jnod).*fun(Inod)).*detJw;    % beta derivative
                
                
                Kxh(:,Inod,Jnod)=Kxh(:,Inod,Jnod)...
                    +(etaint.*(4*exx(:,Iint)+2*eyy(:,Iint)).*Deriv(:,1,Inod).*fun(Jnod)...
                    +etaint.*2.*exy(:,Iint).*Deriv(:,2,Inod).*fun(Jnod)...
                    +deltaint.*beta2int.*uint.*fun(Inod).*fun(Jnod)...
                    +rhog*(Heint+hint.*deltaint).*((rho.*dhdx/rhow-dHdx)*ca-sa).*fun(Inod).*fun(Jnod)...
                    +(rho*rhog/rhow).*hint.*Heint.*ca.*fun(Inod).*Deriv(:,1,Jnod)...
                    -varrhog.*hint.*ca.*Deriv(:,1,Inod).*fun(Jnod)).*detJw;
                
                Kyh(:,Inod,Jnod)=Kyh(:,Inod,Jnod)...
                    +(etaint.*(4*eyy(:,Iint)+2*exx(:,Iint)).*Deriv(:,2,Inod).*fun(Jnod)...
                    +etaint.*2.*exy(:,Iint).*Deriv(:,1,Inod).*fun(Jnod)...
                    +deltaint.*beta2int.*vint.*fun(Inod).*fun(Jnod)...
                    +rhog.*(Heint+hint.*deltaint).*(rho.*dhdy/rhow-dHdy).*ca.*fun(Inod).*fun(Jnod)...
                    +(rho*rhog/rhow)*hint.*Heint.*ca.*fun(Inod).*Deriv(:,2,Jnod)...
                    -varrhog.*hint.*ca.*Deriv(:,2,Inod).*fun(Jnod)).*detJw;
                
                Khu(:,Inod,Jnod)=Khu(:,Inod,Jnod)...
                    +theta*(dhdx.*fun(Jnod)+hint.*Deriv(:,1,Jnod)).*fun(Inod).*detJw;
                
                Khv(:,Inod,Jnod)=Khv(:,Inod,Jnod)...
                    +theta*(dhdy.*fun(Jnod)+hint.*Deriv(:,2,Jnod)).*fun(Inod).*detJw;
                
                Khh(:,Inod,Jnod)=Khh(:,Inod,Jnod)...
                    +(fun(Jnod)/dt...
                    +theta*(exx(:,Iint).*fun(Jnod)+uint.*Deriv(:,1,Jnod))...
                    +theta*(eyy(:,Iint).*fun(Jnod)+vint.*Deriv(:,2,Jnod))).*fun(Inod).*detJw;
                
                
            end
            
        end
        
        t1=-rhog*Heint.*hint.*((rho*dhdx/rhow-dHdx).*ca-sa).*fun(Inod);
        t2=0.5*ca*varrhog.*hint.^2.*Deriv(:,1,Inod);
        t3=hint.*etaint.*(4*exx(:,Iint)+2*eyy(:,Iint)).*Deriv(:,1,Inod);
        t4=hint.*etaint.*2.*exy(:,Iint).*Deriv(:,2,Inod);
        t5=beta2int.*uint.*fun(Inod);
        
        Tx(:,Inod)=Tx(:,Inod)+(t3+t4+t5).*detJw;
        Fx(:,Inod)=Fx(:,Inod)+(t1+t2).*detJw;
        
        % the t1 term should be: rho g h \p_x ( He(h-h_f) (rho h /rho_w-H)
        % which can also be written as: rho g h \p_x s^{'}
        % with s^{'}=He(h-h_f) ( \rho h /\rho h -H)
        % seems to me that the derivative with respect to He is not included in the current version
        
        t1=-rhog.*Heint.*hint.*(rho*dhdy/rhow-dHdy).*ca.*fun(Inod);
        t2=0.5*ca*varrhog.*hint.^2.*Deriv(:,2,Inod);
        t3=hint.*etaint.*(4*eyy(:,Iint)+2*exx(:,Iint)).*Deriv(:,2,Inod);
        t4=hint.*etaint.*2.*exy(:,Iint).*Deriv(:,1,Inod);
        t5=beta2int.*vint.*fun(Inod);
        
        Ty(:,Inod)=Ty(:,Inod)+(t3+t4+t5).*detJw;
        
        Fy(:,Inod)=Fy(:,Inod)+(t1+t2).*detJw;
        
        th=  (theta*(exx(:,Iint).*hint+uint.*dhdx)+(1-theta)*(du0dx.*h0int+u0int.*dh0dx)).*fun(Inod)...
            +(theta*(eyy(:,Iint).*hint+vint.*dhdy)+(1-theta)*(dv0dy.*h0int+v0int.*dh0dy)).*fun(Inod);
        fh=-((hint-h0int)/dt-asint-abint).*fun(Inod);
        
        
        Th(:,Inod)=Th(:,Inod)+th.*detJw;
        Fh(:,Inod)=Fh(:,Inod)+fh.*detJw;
        
        
    end
end



%% assemble right-hand side

T=sparse(neq,1); F=sparse(neq,1);

for Inod=1:nod
    
    
    T=T+sparse(connectivity(:,Inod),ones(Nele,1),Tx(:,Inod),neq,1);
    T=T+sparse(connectivity(:,Inod)+neqx,ones(Nele,1),Ty(:,Inod),neq,1);
    T=T+sparse(connectivity(:,Inod)+2*neqx,ones(Nele,1),Th(:,Inod),neq,1);
    
    F=F+sparse(connectivity(:,Inod),ones(Nele,1),Fx(:,Inod),neq,1);
    F=F+sparse(connectivity(:,Inod)+neqx,ones(Nele,1),Fy(:,Inod),neq,1);
    F=F+sparse(connectivity(:,Inod)+2*neqx,ones(Nele,1),Fh(:,Inod),neq,1);
end
%%

R=T-F;

if ~Ronly
    
    
    
    % large memory version
    largeMemory=1;
    if largeMemory==1
        
        Iind=zeros(9*nod*nod*Nele,1); Jind=zeros(9*nod*nod*Nele,1);Xval=zeros(9*nod*nod*Nele,1);
        istak=0;
        for Inod=1:nod
            for Jnod=1:nod
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod); Jind(istak+1:istak+Nele)=connectivity(:,Jnod); Xval(istak+1:istak+Nele)=Kxu(:,Inod,Jnod);
                istak=istak+Nele;
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod); Jind(istak+1:istak+Nele)=connectivity(:,Jnod)+neqx; Xval(istak+1:istak+Nele)=Kxv(:,Inod,Jnod);
                istak=istak+Nele;
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod); Jind(istak+1:istak+Nele)=connectivity(:,Jnod)+2*neqx; Xval(istak+1:istak+Nele)=Kxh(:,Inod,Jnod);
                istak=istak+Nele;
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod)+neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod); Xval(istak+1:istak+Nele)=Kyu(:,Inod,Jnod);
                istak=istak+Nele;
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod)+neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod)+neqx; Xval(istak+1:istak+Nele)=Kyv(:,Inod,Jnod);
                istak=istak+Nele;
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod)+neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod)+2*neqx; Xval(istak+1:istak+Nele)=Kyh(:,Inod,Jnod);
                istak=istak+Nele;
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod)+2*neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod); Xval(istak+1:istak+Nele)=Khu(:,Inod,Jnod);
                istak=istak+Nele;
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod)+2*neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod)+neqx; Xval(istak+1:istak+Nele)=Khv(:,Inod,Jnod);
                istak=istak+Nele;
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod)+2*neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod)+2*neqx; Xval(istak+1:istak+Nele)=Khh(:,Inod,Jnod);
                istak=istak+Nele;
            end
            
        end
        
        %%
        
        K=sparse(Iind,Jind,Xval,neq,neq);
        
        
        
        
        %nz=nnz(K);
        
        
        %%
        
    else
        Iind=zeros(9*nod*Nele,1); Jind=zeros(9*nod*Nele,1);Xval=zeros(9*nod*Nele,1);
        K=sparse(neq,neq);
        
        for Inod=1:nod
            istak=0;
            for Jnod=1:nod
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod); Jind(istak+1:istak+Nele)=connectivity(:,Jnod); Xval(istak+1:istak+Nele)=Kxu(:,Inod,Jnod);
                istak=istak+Nele;
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod); Jind(istak+1:istak+Nele)=connectivity(:,Jnod)+neqx; Xval(istak+1:istak+Nele)=Kxv(:,Inod,Jnod);
                istak=istak+Nele;
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod); Jind(istak+1:istak+Nele)=connectivity(:,Jnod)+2*neqx; Xval(istak+1:istak+Nele)=Kxh(:,Inod,Jnod);
                istak=istak+Nele;
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod)+neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod); Xval(istak+1:istak+Nele)=Kyu(:,Inod,Jnod);
                istak=istak+Nele;
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod)+neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod)+neqx; Xval(istak+1:istak+Nele)=Kyv(:,Inod,Jnod);
                istak=istak+Nele;
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod)+neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod)+2*neqx; Xval(istak+1:istak+Nele)=Kyh(:,Inod,Jnod);
                istak=istak+Nele;
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod)+2*neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod); Xval(istak+1:istak+Nele)=Khu(:,Inod,Jnod);
                istak=istak+Nele;
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod)+2*neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod)+neqx; Xval(istak+1:istak+Nele)=Khv(:,Inod,Jnod);
                istak=istak+Nele;
                
                Iind(istak+1:istak+Nele)=connectivity(:,Inod)+2*neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod)+2*neqx; Xval(istak+1:istak+Nele)=Khh(:,Inod,Jnod);
                istak=istak+Nele;
            end
            K=K+sparse(Iind,Jind,Xval,neq,neq);
        end
    end
end



end



