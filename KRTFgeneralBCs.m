function [Ruv,Kuv,Tuv,Fuv]=KRTFgeneralBCs(CtrlVar,MUA,s,S,B,h,ub,vb,AGlen,n,C,m,alpha,rho,rhow,g)
                 
% KRTFgeneralBCs(s,S,B,h,u,v,AGlen,n,C,m,coordinates,connectivity,Boundary,nip,alpha,rho,rhow,g,CtrlVar)

if nargout==1
    Ronly=1;
else
    Ronly=0;
end

% calculates the tangent matrix (K) and right-hand side (-R) in a vectorized form

%	Ruv=Tuv-Fuv;
%   Fuv are the `external forces', i.e. right-hand side of the original system
%   Tuv are the `internal forces'. The equation is considered solved once internal and external forces are
%   equal to within a given tolerance

if any(h<0) ; warning('MATLAB:KRTF:hnegative',' h negative ') ; end
if any(C<0) ; warning('MATLAB:KRTF:Cnegative',' C negative ') ; end
if ~isreal(C) ; save TestSave ; error('KRTF: C not real ') ; end

if any(isnan(ub)) ; save TestSave ; error('KRTF: u is nan ') ; end
if any(isnan(vb)) ; save TestSave ; error('KRTF: v is nan ') ; end

if CtrlVar.Piccard
    Dvisk=0;
    Dbeta=0;
else
    Dvisk=CtrlVar.NRviscosity ; % if gradients with respect to visk not to be included set to 0, othewise 1
    Dbeta=CtrlVar.NRbeta2;
end



[etaInt,~,~,exx,eyy,exy,Eint]=calcStrainRatesEtaInt(CtrlVar,MUA,ub,vb,AGlen,n);
if ~isreal(etaInt) ; save TestSave ; error('KRTF: etaInt not real ') ; end
if ~isreal(Eint) ; save TestSave ; error('KRTF: Eint not real ') ; end


%Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
ndim=2; dof=2; neq=dof*MUA.Nnodes;
neqx=MUA.Nnodes ;

%[b,s]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar);

hnod=reshape(h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
snod=reshape(s(MUA.connectivity,1),MUA.Nele,MUA.nod);
unod=reshape(ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(vb(MUA.connectivity,1),MUA.Nele,MUA.nod);


if ~CtrlVar.CisElementBased 
    Cnod=reshape(C(MUA.connectivity,1),MUA.Nele,MUA.nod); 
    mnod=reshape(m(MUA.connectivity,1),MUA.Nele,MUA.nod); 
end


Snod=reshape(S(MUA.connectivity,1),MUA.Nele,MUA.nod);
Bnod=reshape(B(MUA.connectivity,1),MUA.Nele,MUA.nod);
rhonod=reshape(rho(MUA.connectivity,1),MUA.Nele,MUA.nod);


ca=cos(alpha); sa=sin(alpha);


%[points,weights]=sample('triangle',MUA.nip,ndim);

if ~Ronly
    d1d1=zeros(MUA.Nele,MUA.nod,MUA.nod); d2d2=zeros(MUA.Nele,MUA.nod,MUA.nod);  d1d2=zeros(MUA.Nele,MUA.nod,MUA.nod); d2d1=zeros(MUA.Nele,MUA.nod,MUA.nod);
end

Tx=zeros(MUA.Nele,MUA.nod);  Ty=zeros(MUA.Nele,MUA.nod); Fx=zeros(MUA.Nele,MUA.nod);  Fy=zeros(MUA.Nele,MUA.nod);


for Iint=1:MUA.nip
    
    
    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    %[Deriv,detJ]=derivVector(coordinates,connectivity,nip,Iint);
    
    if isfield(MUA,'Deriv') && isfield(MUA,'DetJ') && ~isempty(MUA.Deriv) && ~isempty(MUA.DetJ)
        Deriv=MUA.Deriv(:,:,:,Iint);
        detJ=MUA.DetJ(:,Iint);
    else
        [Deriv,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
    end
    %        fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    %       [Deriv,detJ]=derivVector(coordinates,connectivity,nip,Iint);
    
    % Deriv : Nele x dof x nod
    %  detJ : Nele
    
    % values at integration this point
    hint=hnod*fun;
    sint=snod*fun;
    uint=unod*fun;
    vint=vnod*fun;
    
    if CtrlVar.CisElementBased
        Cint=C;
        mint=m;
    else
        Cint=Cnod*fun;
        Cint(Cint<CtrlVar.Cmin)=CtrlVar.Cmin; % for higher order elements it is possible that Cint is less than any of the nodal values
        mint=mnod*fun;
    end
    
    
    Bint=Bnod*fun;
    Sint=Snod*fun;
    bint=sint-hint;
    Hint=Sint-Bint;
    rhoint=rhonod*fun;
    dint = HeavisideApprox(CtrlVar.kH,Hint,CtrlVar.Hh0).*(Sint-bint);  % draft
    
    hfint=rhow*Hint./rhoint;
    
    Heint = HeavisideApprox(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);
    
    [beta2int,Dbeta2Duuint,Dbeta2Dvvint,Dbeta2Duvint] = calcBeta2in2Dint(uint,vint,Cint,mint,Heint,CtrlVar);
    if ~isreal(beta2int)     ; save KRTFgeneralBCsErrorFile  ; error('KRTF: beta2int not real. All variables saved to ''KRTFgeneralBCsErrorFile'' ') ; end
    if ~isreal(Dbeta2Duuint) ; save KRTFgeneralBCsErrorFile  ; error('KRTF: Dbeta2Duuint not real. All variables saved to ''KRTFgeneralBCsErrorFile'' ') ; end
    if ~isreal(Dbeta2Dvvint) ; save KRTFgeneralBCsErrorFile  ; error('KRTF: Dbeta2Dvvint not real. All variables saved to ''KRTFgeneralBCsErrorFile'' ') ; end
    if ~isreal(Dbeta2Duvint) ; save KRTFgeneralBCsErrorFile  ; error('KRTF: Dbeta2Duvint not real. All variables saved to ''KRTFgeneralBCsErrorFile'' ') ; end
    
    
    etaint=etaInt(:,Iint) ;  % I could consider calculating this here
    
    
    % derivatives at this integration point for all elements
    dsdx=zeros(MUA.Nele,1); dhdx=zeros(MUA.Nele,1);
    dsdy=zeros(MUA.Nele,1); dhdy=zeros(MUA.Nele,1);
    
    
    
    
    for Inod=1:MUA.nod
        
        dsdx=dsdx+Deriv(:,1,Inod).*snod(:,Inod);
        dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
        dsdy=dsdy+Deriv(:,2,Inod).*snod(:,Inod);
        dhdy=dhdy+Deriv(:,2,Inod).*hnod(:,Inod);
        
    end
    
    dbdx=dsdx-dhdx; dbdy=dsdy-dhdy;
    
    detJw=detJ*MUA.weights(Iint);
    
    
    for Inod=1:MUA.nod
        if ~Ronly
            for Jnod=1:MUA.nod
                
                
                d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)...
                    +(4*hint.*etaint.*Deriv(:,1,Inod).*Deriv(:,1,Jnod)...
                    +hint.*etaint.*Deriv(:,2,Inod).*Deriv(:,2,Jnod)...
                    +beta2int           .*fun(Jnod).*fun(Inod)...    % basal friction, Weertman
                    +Dbeta.*Dbeta2Duuint.*fun(Jnod).*fun(Inod))...   % basal friction, Weertman, directional derivative, uu
                    .*detJw;  
                
                
                d2d2(:,Inod,Jnod)=d2d2(:,Inod,Jnod)...
                    +(4*hint.*etaint.*Deriv(:,2,Inod).*Deriv(:,2,Jnod)...
                    +hint.*etaint.*Deriv(:,1,Inod).*Deriv(:,1,Jnod)...
                    +beta2int           .*fun(Jnod).*fun(Inod)...  % basal friction, Weertman
                    +Dbeta.*Dbeta2Dvvint.*fun(Jnod).*fun(Inod))... % basal friction, Weertman, directional derivative, vv
                    .*detJw ;
               
                
                
                %+Dbeta*vint.*Dbeta2Dvint.*fun(Jnod).*fun(Inod)).*detJw ;
                
                
                d1d2(:,Inod,Jnod)=d1d2(:,Inod,Jnod)...
                    +(etaint.*hint.*(2*Deriv(:,1,Inod).*Deriv(:,2,Jnod)+Deriv(:,2,Inod).*Deriv(:,1,Jnod))...
                    +Dbeta.*Dbeta2Duvint.*fun(Jnod).*fun(Inod))...    % beta derivative, uv
                    .*detJw; 
                
                
                d2d1(:,Inod,Jnod)=d2d1(:,Inod,Jnod)...
                    +(etaint.*hint.*(2*Deriv(:,2,Inod).*Deriv(:,1,Jnod)+Deriv(:,1,Inod).*Deriv(:,2,Jnod))...
                    +Dbeta.*Dbeta2Duvint*fun(Jnod).*fun(Inod)).*detJw;    % beta derivative, uv
                
                
                %                dxu=E (2 exx+eyy)
                %                dyu=E exy
                %                dyv=E (2 eyy + exx )
                %                dxv=E exy = dyu
                
                % 				xy12J=(2*exx(:,Iint)+eyy(:,Iint)).*Deriv(:,1,Jnod)+exy(:,Iint).*Deriv(:,2,Jnod);
                % 				yx21J=(2*eyy(:,Iint)+exx(:,Iint)).*Deriv(:,2,Jnod)+exy(:,Iint).*Deriv(:,1,Jnod);
                %
                % 				xy12I=(2*exx(:,Iint)+eyy(:,Iint)).*Deriv(:,1,Inod)+exy(:,Iint).*Deriv(:,2,Inod);
                % 				yx21I=(2*eyy(:,Iint)+exx(:,Iint)).*Deriv(:,2,Inod)+exy(:,Iint).*Deriv(:,1,Inod);
                %
                %
                %
                % 				E11=2*Eint(:,Iint).*hint.*xy12J.*xy12I;
                % 				E12=2*Eint(:,Iint).*hint.*yx21J.*xy12I;
                % 				E21=2*Eint(:,Iint).*hint.*xy12J.*yx21I;
                % 				E22=2*Eint(:,Iint).*hint.*yx21I.*yx21J;
                
                Deu=Eint(:,Iint).*((2*exx(:,Iint)+eyy(:,Iint)).*Deriv(:,1,Jnod)+exy(:,Iint).*Deriv(:,2,Jnod));
                Dev=Eint(:,Iint).*((2*eyy(:,Iint)+exx(:,Iint)).*Deriv(:,2,Jnod)+exy(:,Iint).*Deriv(:,1,Jnod));
                
                % E11=h Deu (4 p_x u + 2 p_y v)   + h Deu  ( p_x v + p_y u) p_y N_p
                
                E11=  hint.*(4.*exx(:,Iint)+2.*eyy(:,Iint)).*Deu.*Deriv(:,1,Inod)...
                    +2*hint.*exy(:,Iint).*Deu.*Deriv(:,2,Inod);
                
                
                E12=  hint.*(4.*exx(:,Iint)+2.*eyy(:,Iint)).*Dev.*Deriv(:,1,Inod)...
                    +2*hint.*exy(:,Iint).*Dev.*Deriv(:,2,Inod);
                
                
                
                E22=  hint.*(4.*eyy(:,Iint)+2.*exx(:,Iint)).*Dev.*Deriv(:,2,Inod)...
                    +2*hint.*exy(:,Iint).*Dev.*Deriv(:,1,Inod);
                
                
                E21= hint.*(4.*eyy(:,Iint)+2.*exx(:,Iint)).*Deu.*Deriv(:,2,Inod)...
                    +2*hint.*exy(:,Iint).*Deu.*Deriv(:,1,Inod);
                
                
                
                d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+Dvisk*E11.*detJw;
                d2d2(:,Inod,Jnod)=d2d2(:,Inod,Jnod)+Dvisk*E22.*detJw;
                d1d2(:,Inod,Jnod)=d1d2(:,Inod,Jnod)+Dvisk*E12.*detJw;
                d2d1(:,Inod,Jnod)=d2d1(:,Inod,Jnod)+Dvisk*E21.*detJw;
                
            end
        end
        
        t1=-g*(rhoint.*hint-rhow*dint).*dbdx.*fun(Inod)*ca+ rhoint.*g.*hint.*sa.*fun(Inod);
        
        t2=0.5*ca*g.*(rhoint.*hint.^2-rhow.*dint.^2).*Deriv(:,1,Inod);
        t3=hint.*etaint.*(4*exx(:,Iint)+2*eyy(:,Iint)).*Deriv(:,1,Inod);
        t4=hint.*etaint.*2.*exy(:,Iint).*Deriv(:,2,Inod);
        t5=beta2int.*uint.*fun(Inod);  % basal friction, Weertman, u
        
        Tx(:,Inod)=Tx(:,Inod)+(t3+t4+t5).*detJw;
        Fx(:,Inod)=Fx(:,Inod)+(t1+t2).*detJw;
        
        t1=-g*(rhoint.*hint-rhow*dint).*dbdy.*fun(Inod)*ca;
        t2=0.5*ca*g.*(rhoint.*hint.^2-rhow.*dint.^2).*Deriv(:,2,Inod);
        t3=hint.*etaint.*(4*eyy(:,Iint)+2*exx(:,Iint)).*Deriv(:,2,Inod);
        t4=hint.*etaint.*2.*exy(:,Iint).*Deriv(:,1,Inod);
        t5=beta2int.*vint.*fun(Inod); % basal friction, Weertman, v
        
        Ty(:,Inod)=Ty(:,Inod)+(t3+t4+t5).*detJw;
        Fy(:,Inod)=Fy(:,Inod)+(t1+t2).*detJw;
        
        
        
        
        
    end
end

% add boundary integral relatedto Dirichlet boundary conditions


% assemble right-hand side

Tuv=sparseUA(neq,1); Fuv=sparseUA(neq,1);

for Inod=1:MUA.nod
    
    
    Tuv=Tuv+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),Tx(:,Inod),neq,1);
    Tuv=Tuv+sparseUA(MUA.connectivity(:,Inod)+neqx,ones(MUA.Nele,1),Ty(:,Inod),neq,1);
    
    Fuv=Fuv+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),Fx(:,Inod),neq,1);
    Fuv=Fuv+sparseUA(MUA.connectivity(:,Inod)+neqx,ones(MUA.Nele,1),Fy(:,Inod),neq,1);
end

Ruv=Tuv-Fuv;

if ~Ronly
    iSparse=1;  %	faster
    if iSparse==1
        % uses the sparse function less often
        
        
        %Iind=zeros(nod*Nele*4,1); Jind=zeros(nod*Nele*4,1);Xval=zeros(nod*Nele*4,1);
        Iind=zeros(MUA.nod*MUA.nod*MUA.Nele*4,1); Jind=zeros(MUA.nod*MUA.nod*MUA.Nele*4,1);Xval=zeros(MUA.nod*MUA.nod*MUA.Nele*4,1); istak=0;
        for Inod=1:MUA.nod
            %istak=0;
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
            %K=K+sparse(Iind,Jind,Xval,neq,neq);
        end
        % nzmax=size(unique([Iind Jind],'rows'),1) ; K=sparse(Iind,Jind,Xval,neq,neq,nzmax); not sure why this
        % does not work
        
        Kuv=sparseUA(Iind,Jind,Xval,neq,neq);
        
    else
        % creates the sparse matrix in steps, requires no extra
        for Inod=1:MUA.nod
            for Jnod=1:MUA.nod
                Kuv=Kuv+sparseUA(MUA.connectivity(:,Inod),MUA.connectivity(:,Jnod),d1d1(:,Inod,Jnod),neq,neq);
                Kuv=Kuv+sparseUA(MUA.connectivity(:,Inod)+neqx,MUA.connectivity(:,Jnod)+neqx,d2d2(:,Inod,Jnod),neq,neq);
                Kuv=Kuv+sparseUA(MUA.connectivity(:,Inod),MUA.connectivity(:,Jnod)+neqx,d1d2(:,Inod,Jnod),neq,neq);
                Kuv=Kuv+sparseUA(MUA.connectivity(:,Inod)+neqx,MUA.connectivity(:,Jnod),d2d1(:,Inod,Jnod),neq,neq);
                %K=K+sparse(connectivity(:,Inod)+neqx,connectivity(:,Jnod),d1d2(:,Jnod,Inod),neq,neq);
            end
        end
    end
    
    
    Kuv=(Kuv+Kuv.')/2 ; % I know that the matrix must be symmetric, but numerically this may not be strickly so
    % Note: for numerical verification of distributed parameter gradient it is important to
    % not to use the complex conjugate transpose.
    % whos('K')
    
    % Boundary contribution
    
    if CtrlVar.IncludeDirichletBoundaryIntegralDiagnostic
        [KBoundary,rhsBoundary]=DirichletBoundaryIntegralDiagnostic(MUA.coordinates,MUA.connectivity,Boundary,nip,h,ub,vb,AGlen,n,alpha,rho,rhow,g,CtrlVar);
        Kuv=Kuv+KBoundary ; Ruv=Ruv+rhsBoundary;
    end
end
end



