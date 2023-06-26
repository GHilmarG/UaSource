function [Ruv,Kuv,Tint,Fext]=uvMatrixAssemblySSTREAM(CtrlVar,MUA,F)

%
% Ruv=Tint-Fext;
% Tint   : internal nodal forces
% Fint   : external nodal forces

narginchk(3,5)
nargoutchk(1,5)



ZeroFields=CtrlVar.uvAssembly.ZeroFields;
Ronly=CtrlVar.uvMatrixAssembly.Ronly;

if Ronly
    Kuv=[];
end



if ZeroFields
    F.ub=F.ub*0;
    F.vb=F.vb*0;
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


if CtrlVar.Picard
    Dvisk=0;
    Dbeta=0;
else
    Dvisk=CtrlVar.NRviscosity ; % if gradients with respect to visk not to be included set to 0, otherwise 1
    Dbeta=CtrlVar.NRbeta2;
end

g=F.g;





%Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
ndim=2; dof=2; neq=dof*MUA.Nnodes;
neqx=MUA.Nnodes ;

%[b,s]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar);

hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
snod=reshape(F.s(MUA.connectivity,1),MUA.Nele,MUA.nod);

ubnod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vbnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);


H=F.S-F.B;


if CtrlVar.IncludeMelangeModelPhysics
    
    uonod=reshape(F.uo(MUA.connectivity,1),MUA.Nele,MUA.nod);
    vonod=reshape(F.vo(MUA.connectivity,1),MUA.Nele,MUA.nod);
    
    uanod=reshape(F.ua(MUA.connectivity,1),MUA.Nele,MUA.nod);
    vanod=reshape(F.va(MUA.connectivity,1),MUA.Nele,MUA.nod);

end

%if ~CtrlVar.CisElementBased

Cnod=reshape(F.C(MUA.connectivity,1),MUA.Nele,MUA.nod);
mnod=reshape(F.m(MUA.connectivity,1),MUA.Nele,MUA.nod);

if ~isempty(F.q)
    qnod=reshape(F.q(MUA.connectivity,1),MUA.Nele,MUA.nod);
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


%if ~CtrlVar.AGlenisElementBased
AGlennod=reshape(F.AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
nnod=reshape(F.n(MUA.connectivity,1),MUA.Nele,MUA.nod);
%end



Snod=reshape(F.S(MUA.connectivity,1),MUA.Nele,MUA.nod);
Bnod=reshape(F.B(MUA.connectivity,1),MUA.Nele,MUA.nod);
Hnod=Snod-Bnod;
rhonod=reshape(F.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);

ca=cos(F.alpha); sa=sin(F.alpha);


if CtrlVar.uvGroupAssembly
    hfnod=F.rhow*(Snod-Bnod)./rhonod;
    bnod=reshape(F.b(MUA.connectivity,1),MUA.Nele,MUA.nod);
    dnod = HeavisideApprox(CtrlVar.kH,Hnod,CtrlVar.Hh0).*(Snod-bnod);  % draft
    deltanod=DiracDelta(CtrlVar.kH,hnod-hfnod,CtrlVar.Hh0);
    Henod = HeavisideApprox(CtrlVar.kH,hnod-hfnod,CtrlVar.Hh0);
end


%[points,weights]=sample('triangle',MUA.nip,ndim);

if ~Ronly
    d1d1=zeros(MUA.Nele,MUA.nod,MUA.nod); d2d2=zeros(MUA.Nele,MUA.nod,MUA.nod);  d1d2=zeros(MUA.Nele,MUA.nod,MUA.nod); d2d1=zeros(MUA.Nele,MUA.nod,MUA.nod);
end

Tx=zeros(MUA.Nele,MUA.nod);  Ty=zeros(MUA.Nele,MUA.nod); Fx=zeros(MUA.Nele,MUA.nod);  Fy=zeros(MUA.Nele,MUA.nod);


for Iint=1:MUA.nip


    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points

    % if isfield(MUA,'Deriv') && isfield(MUA,'DetJ') && ~isempty(MUA.Deriv) && ~isempty(MUA.DetJ)
    Deriv=MUA.Deriv(:,:,:,Iint);  % Deriv at integration points
    detJ=MUA.DetJ(:,Iint);
    % else
    %     [Deriv,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,MUA.points,Iint);
    % end

    
    %        fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    %       [Deriv,detJ]=derivVector(coordinates,connectivity,nip,Iint);
    
    % Deriv : Nele x dof x nod
    %  detJ : Nele
    
    % values at integration this point
    hint=hnod*fun;
    sint=snod*fun;
    
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




    if CtrlVar.IncludeMelangeModelPhysics
        Coint=Conod*fun;
        moint=monod*fun;

        Caint=Canod*fun;
        maint=manod*fun;
    end
    %   end


    % if CtrlVar.AGlenisElementBased
    %     AGlenint=F.AGlen;
    %     nint=F.n;
    % else
        AGlenint=AGlennod*fun;
        AGlenint(AGlenint<CtrlVar.AGlenmin)=CtrlVar.AGlenmin;
        nint=nnod*fun;
   %  end



    
    Bint=Bnod*fun;
    Sint=Snod*fun;
    bint=sint-hint;
    Hint=Sint-Bint;
    rhoint=rhonod*fun;
    

  
    %
    

    % deltaint=DiracDelta(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);      
    % Heint = HeavisideApprox(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);
    
    

    if CtrlVar.uvGroupAssembly
        %% interpolating dint, hfint, Heint and deltaint onto the    integration points
        dint=dnod*fun;
        deltaint=deltanod*fun;
        Heint=Henod*fun;

      % if any(dint<0)
      %     fprintf(" dint negative \n")
      % end
      % if any(deltaint<0)
      %     fprintf(" deltaint negative \n")
      % end
      % if any(Heint<0)
      %     fprintf(" Heint negative \n")
      % end

    else
        %% evaluating dint, hfint, Heint and deltaint at integration points#
        hfint=F.rhow*Hint./rhoint;
        dint = HeavisideApprox(CtrlVar.kH,Hint,CtrlVar.Hh0).*(Sint-bint);  % draft
        
        
        Heint = HeavisideApprox(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);
        deltaint=DiracDelta(CtrlVar.kH,hint-hfint,CtrlVar.Hh0); % dHeint/dh
        
    end
    
    % derivatives at this integration point for all elements
    dsdx=zeros(MUA.Nele,1); dhdx=zeros(MUA.Nele,1);
    dsdy=zeros(MUA.Nele,1); dhdy=zeros(MUA.Nele,1);
    
    exx=zeros(MUA.Nele,1);
    eyy=zeros(MUA.Nele,1);
    exy=zeros(MUA.Nele,1);
    
    
    for Inod=1:MUA.nod
        
        dsdx=dsdx+Deriv(:,1,Inod).*snod(:,Inod);
        dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
        dsdy=dsdy+Deriv(:,2,Inod).*snod(:,Inod);
        dhdy=dhdy+Deriv(:,2,Inod).*hnod(:,Inod);
        
        exx=exx+Deriv(:,1,Inod).*ubnod(:,Inod);
        eyy=eyy+Deriv(:,2,Inod).*vbnod(:,Inod);
        exy=exy+0.5*(Deriv(:,1,Inod).*vbnod(:,Inod) + Deriv(:,2,Inod).*ubnod(:,Inod));
        
        
    end
    
    

    [taux,tauy,dtauxdu,dtauxdv,dtauydu,dtauydv] = ...
        BasalDrag(CtrlVar,MUA,Heint,deltaint,hint,Bint,Hint,rhoint,F.rhow,uint,vint,Cint,mint,uoint,voint,Coint,moint,uaint,vaint,Caint,maint,qint,g,mukint);
    [etaint,Eint]=EffectiveViscositySSTREAM(CtrlVar,AGlenint,nint,exx,eyy,exy);

  

    dbdx=dsdx-dhdx; dbdy=dsdy-dhdy;
    
    detJw=detJ*MUA.weights(Iint);
    
    
    for Inod=1:MUA.nod
        if ~Ronly
            for Jnod=1:MUA.nod
                
                
                d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)...
                    +(4*hint.*etaint.*Deriv(:,1,Inod).*Deriv(:,1,Jnod)...
                    +hint.*etaint.*Deriv(:,2,Inod).*Deriv(:,2,Jnod)...
                    +dtauxdu.*fun(Jnod).*fun(Inod)... %+beta2int.*fun(Jnod).*fun(Inod)+Dbeta.*Dbeta2Duuint.*fun(Jnod).*fun(Inod))...   
                    ).*detJw;  
                
                
                d2d2(:,Inod,Jnod)=d2d2(:,Inod,Jnod)...
                    +(4*hint.*etaint.*Deriv(:,2,Inod).*Deriv(:,2,Jnod)...
                    +hint.*etaint.*Deriv(:,1,Inod).*Deriv(:,1,Jnod)...
                    +dtauydv.*fun(Jnod).*fun(Inod)...   %+beta2int.*fun(Jnod).*fun(Inod)+Dbeta.*Dbeta2Dvvint.*fun(Jnod).*fun(Inod))... 
                    ).*detJw ;
               
                
                
                d1d2(:,Inod,Jnod)=d1d2(:,Inod,Jnod)...
                    +(etaint.*hint.*(2*Deriv(:,1,Inod).*Deriv(:,2,Jnod)+Deriv(:,2,Inod).*Deriv(:,1,Jnod))...
                    + +dtauxdv.*fun(Jnod).*fun(Inod)...   % Dbeta.*Dbeta2Duvint.*fun(Jnod).*fun(Inod))...    % beta derivative, uv
                    ).*detJw;
                
                
                d2d1(:,Inod,Jnod)=d2d1(:,Inod,Jnod)...
                    +(etaint.*hint.*(2*Deriv(:,2,Inod).*Deriv(:,1,Jnod)+Deriv(:,1,Inod).*Deriv(:,2,Jnod))...
                    +dtauydu.*fun(Jnod).*fun(Inod)...    %+Dbeta.*Dbeta2Duvint*fun(Jnod).*fun(Inod)).*detJw;    % beta derivative, uv
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
                
                
                
                d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+Dvisk*E11.*detJw;
                d2d2(:,Inod,Jnod)=d2d2(:,Inod,Jnod)+Dvisk*E22.*detJw;
                d1d2(:,Inod,Jnod)=d1d2(:,Inod,Jnod)+Dvisk*E12.*detJw;
                d2d1(:,Inod,Jnod)=d2d1(:,Inod,Jnod)+Dvisk*E21.*detJw;
                
            end
        end
        
        t1=-F.g*(rhoint.*hint-F.rhow*dint).*dbdx.*fun(Inod)*ca+ rhoint.*F.g.*hint.*sa.*fun(Inod);
        t2=0.5*F.g.*ca*(rhoint.*hint.^2-F.rhow.*dint.^2).*Deriv(:,1,Inod);

        % 0.5*F.g.*ca* rhoint.*hint.^2 .*Deriv(:,1,Inod);  +   0.5*F.g.*ca* rhoint.*hint.^2 .*fun(Inod).*drhodx
        
        t3=hint.*etaint.*(4*exx+2*eyy).*Deriv(:,1,Inod);
        t4=hint.*etaint.*2.*exy.*Deriv(:,2,Inod);
        t5=taux.*fun(Inod); 
        
        Tx(:,Inod)=Tx(:,Inod)+(t3+t4+t5).*detJw;
        Fx(:,Inod)=Fx(:,Inod)+(t1+t2).*detJw;
        
        t1=-F.g*ca*(rhoint.*hint-F.rhow*dint).*dbdy.*fun(Inod);
        t2=0.5*ca*F.g.*(rhoint.*hint.^2-F.rhow.*dint.^2).*Deriv(:,2,Inod);
        
        t3=hint.*etaint.*(4*eyy+2*exx).*Deriv(:,2,Inod);
        t4=hint.*etaint.*2.*exy.*Deriv(:,1,Inod);
        t5=tauy.*fun(Inod); 
        
        Ty(:,Inod)=Ty(:,Inod)+(t3+t4+t5).*detJw;
        Fy(:,Inod)=Fy(:,Inod)+(t1+t2).*detJw;
        
        
        
        
        
    end
end

% add boundary integral relatedto Dirichlet boundary conditions


% assemble right-hand side

Tint=sparseUA(neq,1); Fext=sparseUA(neq,1);

for Inod=1:MUA.nod
    
    
    Tint=Tint+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),Tx(:,Inod),neq,1);
    Tint=Tint+sparseUA(MUA.connectivity(:,Inod)+neqx,ones(MUA.Nele,1),Ty(:,Inod),neq,1);
    
    Fext=Fext+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),Fx(:,Inod),neq,1);
    Fext=Fext+sparseUA(MUA.connectivity(:,Inod)+neqx,ones(MUA.Nele,1),Fy(:,Inod),neq,1);
end

Ruv=Tint-Fext;

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
    
    if CtrlVar.TestForRealValues
        Kuv=(Kuv+Kuv.')/2 ;
    else
        Kuv=(Kuv+Kuv')/2 ;
    end
    
    
    % I know that the matrix must be symmetric, but numerically this may not be strickly so
    % Note: for numerical verification of distributed parameter gradient it is important to
    % not to use the complex conjugate transpose.
    % whos('K')
    
    % Boundary contribution
    
%     if CtrlVar.IncludeDirichletBoundaryIntegralDiagnostic
%         [KBoundary,rhsBoundary]=DirichletBoundaryIntegralDiagnostic(MUA.coordinates,MUA.connectivity,Boundary,nip,h,ub,vb,AGlen,n,alpha,rho,rhow,g,CtrlVar);
%         Kuv=Kuv+KBoundary ; Ruv=Ruv+rhsBoundary;
%     end
end




end



