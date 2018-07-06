function [Kuv,Ruv,Tuv,Fuv]=KRTF(s,S,B,h,u,v,AGlen,n,C,m,coordinates,connectivity,Boundary,nip,alpha,rho,rhow,g,CtrlVar)
	               
	
	% calculates the tangent matrix (K) and right-hand side (-R) in a vectorized form

	%	R=T-F;
	%   F are the `external forces', i.e. right-hand side of the original system
	%   T are the `internal forces'. The equation is considered solved once internal and external forces are
	%   equal to within a given tolerance 
		
    if any(h<0) ; warning('MATLAB:KRTF:hnegative',' h negative ') ; end
	if any(C<0) ; warning('MATLAB:KRTF:Cnegative',' C negative ') ; end
	if ~isreal(C) ; error('KRTF: C not real ') ; end
	
	
	Dvisk=CtrlVar.NRviscosity ; % if gradients with respect to visk not to be included set to 0, othewise 1
	Dbeta=CtrlVar.NRbeta2;

 
	[etaInt,~,~,exx,eyy,exy,Eint]=calcStrainRatesEtaInt(u,v,coordinates,connectivity,nip,AGlen,n,CtrlVar);
	if ~isreal(etaInt) ; error('KRTF: etaInt not real ') ; end
	if ~isreal(Eint) ; error('KRTF: Eint not real ') ; end
	
	
	Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
	ndim=2; dof=2; neq=dof*Nnodes;
    neqx=Nnodes ;
    
    H=S-B;
    hf=rhow*H./rho ;
    He = HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);  % 1 if grounded, 0 if afloat
    
    %[b,s]=Calc_bs_From_hBS(h,S,B,rho,rhow,CtrlVar);
    
    sdot=He.*(rho.*h/rhow-H); %  sdotTest=s-S-(1-rho/rhow).*h;  norm(sdot-sdotTest) give zero
    
    
    
	hnod=reshape(h(connectivity,1),Nele,nod);   % Nele x nod
	snod=reshape(s(connectivity,1),Nele,nod);
	unod=reshape(u(connectivity,1),Nele,nod);
	vnod=reshape(v(connectivity,1),Nele,nod);
	Cnod=reshape(C(connectivity,1),Nele,nod);
	Snod=reshape(S(connectivity,1),Nele,nod);
    Bnod=reshape(B(connectivity,1),Nele,nod);
	sdotnod=reshape(sdot(connectivity,1),Nele,nod);
    rhonod=reshape(rho(connectivity,1),Nele,nod);
    
	rhog=rho*g; ca=cos(alpha); sa=sin(alpha);
	
	
	[points,weights]=sample('triangle',nip,ndim);

	
	d1d1=zeros(Nele,nod,nod); d2d2=zeros(Nele,nod,nod);  d1d2=zeros(Nele,nod,nod); d2d1=zeros(Nele,nod,nod);
	
	Tx=zeros(Nele,nod);  Ty=zeros(Nele,nod); Fx=zeros(Nele,nod);  Fy=zeros(Nele,nod);
	
    kH=CtrlVar.kH;
    
	for Iint=1:nip
		
				
		fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
		[Deriv,detJ]=derivVector(coordinates,connectivity,nip,Iint);
		% Deriv : Nele x dof x nod
		%  detJ : Nele
		
		% values at integration this point
		hint=hnod*fun;
		%sint=snod*fun;
		uint=unod*fun;
		vint=vnod*fun;
		Cint=Cnod*fun;
		Bint=Bnod*fun; 
        Sint=Snod*fun;
        rhoint=rhonod*fun;
        
        hfint=rhow*(Sint-Bint)./rhoint;
		%Heint=hint > hfint;  % not a good idea
		
		Heint = HeavisideApprox(kH,hint-hfint,CtrlVar.Hh0);  
        
        [beta2int,Dbeta2Duuint,Dbeta2Dvvint,Dbeta2Duvint] = calcBeta2in2Dint(uint,vint,Cint,m,Heint,CtrlVar);
		if ~isreal(beta2int) ; save TestFile beta2int ; error('KRTF: beta2int not real ') ; end
		if ~isreal(Dbeta2Duuint) ; error('KRTF: Dbeta2Duuint not real ') ; end
		if ~isreal(Dbeta2Dvvint) ; error('KRTF: Dbeta2Dvvint not real ') ; end
		if ~isreal(Dbeta2Duvint) ; error('KRTF: Dbeta2Duvint not real ') ; end
		
		
		etaint=etaInt(:,Iint) ;  % I could consider calculating this here
		
		
		% derivatives at this integration point for all elements
		dsdx=zeros(Nele,1); dhdx=zeros(Nele,1); dsdy=zeros(Nele,1); dhdy=zeros(Nele,1);
        dsdotdx=zeros(Nele,1); dsdotdy=zeros(Nele,1);
		
		
		
		for Inod=1:nod
			
            dsdx=dsdx+Deriv(:,1,Inod).*snod(:,Inod);
            dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
            dsdy=dsdy+Deriv(:,2,Inod).*snod(:,Inod);
            dhdy=dhdy+Deriv(:,2,Inod).*hnod(:,Inod);
            dsdotdx=dsdotdx+Deriv(:,1,Inod).*sdotnod(:,Inod);
            dsdotdy=dsdotdy+Deriv(:,2,Inod).*sdotnod(:,Inod);
			
		end

		detJw=detJ*weights(Iint);
		
		
		for Inod=1:nod
			for Jnod=1:nod
				
				
				d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)...
					+(4*hint.*etaint.*Deriv(:,1,Inod).*Deriv(:,1,Jnod)...
					+hint.*etaint.*Deriv(:,2,Inod).*Deriv(:,2,Jnod)...
					+beta2int.*fun(Jnod).*fun(Inod)...
					+Dbeta.*Dbeta2Duuint.*fun(Jnod).*fun(Inod)).*detJw;  % beta derivative	

				
				d2d2(:,Inod,Jnod)=d2d2(:,Inod,Jnod)...
					+(4*hint.*etaint.*Deriv(:,2,Inod).*Deriv(:,2,Jnod)...
					+hint.*etaint.*Deriv(:,1,Inod).*Deriv(:,1,Jnod)...
					+beta2int.*fun(Jnod).*fun(Inod)...
					+Dbeta.*Dbeta2Dvvint.*fun(Jnod).*fun(Inod)).*detJw ;
				%+Dbeta*vint.*Dbeta2Dvint.*fun(Jnod).*fun(Inod)).*detJw ;
				
				
				d1d2(:,Inod,Jnod)=d1d2(:,Inod,Jnod)...
					+(etaint.*hint.*(2*Deriv(:,1,Inod).*Deriv(:,2,Jnod)+Deriv(:,2,Inod).*Deriv(:,1,Jnod))...
					+Dbeta.*Dbeta2Duvint.*fun(Jnod).*fun(Inod)).*detJw;    % beta derivative
				
				
				d2d1(:,Inod,Jnod)=d2d1(:,Inod,Jnod)...
					+(etaint.*hint.*(2*Deriv(:,2,Inod).*Deriv(:,1,Jnod)+Deriv(:,1,Inod).*Deriv(:,2,Jnod))...
					+Dbeta.*Dbeta2Duvint*fun(Jnod).*fun(Inod)).*detJw;    % beta derivative
				
				
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
			
            t1=-g*rhoint.*hint.*(dsdotdx*ca-sa).*fun(Inod);
			%t1=-rhog*hint.*((dsdx-(1-rho/rhow).*dhdx)*ca-sa).*fun(Inod);
			t2=0.5*g*ca*rhoint.*(1-rhoint/rhow).*hint.^2.*Deriv(:,1,Inod);
			t3=hint.*etaint.*(4*exx(:,Iint)+2*eyy(:,Iint)).*Deriv(:,1,Inod);
			t4=hint.*etaint.*2.*exy(:,Iint).*Deriv(:,2,Inod);
			t5=beta2int.*uint.*fun(Inod);
			
			Tx(:,Inod)=Tx(:,Inod)+(t3+t4+t5).*detJw;
			Fx(:,Inod)=Fx(:,Inod)+(t1+t2).*detJw;
			
			
			t1=-g*rhoint.*hint.*dsdotdy*ca.*fun(Inod);
			%t1=-rhog.*hint.*(dsdy-(1-rho/rhow).*dhdy).*ca.*fun(Inod);
			t2=0.5*g*ca*rhoint.*(1-rhoint/rhow).*hint.^2.*Deriv(:,2,Inod);
			t3=hint.*etaint.*(4*eyy(:,Iint)+2*exx(:,Iint)).*Deriv(:,2,Inod);
			t4=hint.*etaint.*2.*exy(:,Iint).*Deriv(:,1,Inod);
			t5=beta2int.*vint.*fun(Inod);
			
			Ty(:,Inod)=Ty(:,Inod)+(t3+t4+t5).*detJw;
			Fy(:,Inod)=Fy(:,Inod)+(t1+t2).*detJw;
			
			
			
			
			
		end
	end
	
	% add boundary integratl relatedto Dirichlet boundary conditions
	
	
	% assemble right-hand side
	
	Tuv=sparse(neq,1); Fuv=sparse(neq,1);
	
	for Inod=1:nod
		
		
		Tuv=Tuv+sparse(connectivity(:,Inod),ones(Nele,1),Tx(:,Inod),neq,1);
		Tuv=Tuv+sparse(connectivity(:,Inod)+neqx,ones(Nele,1),Ty(:,Inod),neq,1);
		
		Fuv=Fuv+sparse(connectivity(:,Inod),ones(Nele,1),Fx(:,Inod),neq,1);
		Fuv=Fuv+sparse(connectivity(:,Inod)+neqx,ones(Nele,1),Fy(:,Inod),neq,1);
	end
	
	Ruv=Tuv-Fuv;
	
	%K=spalloc(neq,neq,neq*25); % not sure how much space is needed but this seems reasonable
	%K=sparse(neq,neq);
	%whos('K')
	% assemble matrix
	iSparse=1;  %	faster
	if iSparse==1
		% uses the sparse function less often
		
		
		%Iind=zeros(nod*Nele*4,1); Jind=zeros(nod*Nele*4,1);Xval=zeros(nod*Nele*4,1);
		Iind=zeros(nod*nod*Nele*4,1); Jind=zeros(nod*nod*Nele*4,1);Xval=zeros(nod*nod*Nele*4,1); istak=0;
		for Inod=1:nod
			%istak=0;
			for Jnod=1:nod
				
				Iind(istak+1:istak+Nele)=connectivity(:,Inod); Jind(istak+1:istak+Nele)=connectivity(:,Jnod); Xval(istak+1:istak+Nele)=d1d1(:,Inod,Jnod);
				istak=istak+Nele;
				
				Iind(istak+1:istak+Nele)=connectivity(:,Inod)+neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod)+neqx; Xval(istak+1:istak+Nele)=d2d2(:,Inod,Jnod);
				istak=istak+Nele;
				
				Iind(istak+1:istak+Nele)=connectivity(:,Inod); Jind(istak+1:istak+Nele)=connectivity(:,Jnod)+neqx; Xval(istak+1:istak+Nele)=d1d2(:,Inod,Jnod);
				istak=istak+Nele;
				
				Iind(istak+1:istak+Nele)=connectivity(:,Inod)+neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod); Xval(istak+1:istak+Nele)=d2d1(:,Inod,Jnod);
				%Iind(istak+1:istak+Nele)=connectivity(:,Inod)+neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod); Xval(istak+1:istak+Nele)=d1d2(:,Jnod,Inod);
				istak=istak+Nele;
				
			end
			%K=K+sparse(Iind,Jind,Xval,neq,neq);
		end
		% nzmax=size(unique([Iind Jind],'rows'),1) ; K=sparse(Iind,Jind,Xval,neq,neq,nzmax); not sure why this
		% does not work
				
		Kuv=sparse(Iind,Jind,Xval,neq,neq);
		
	else
		% creates the sparse matrix in steps, requires no extra
		for Inod=1:nod
			for Jnod=1:nod
				Kuv=Kuv+sparse(connectivity(:,Inod),connectivity(:,Jnod),d1d1(:,Inod,Jnod),neq,neq);
				Kuv=Kuv+sparse(connectivity(:,Inod)+neqx,connectivity(:,Jnod)+neqx,d2d2(:,Inod,Jnod),neq,neq);
				Kuv=Kuv+sparse(connectivity(:,Inod),connectivity(:,Jnod)+neqx,d1d2(:,Inod,Jnod),neq,neq);
				Kuv=Kuv+sparse(connectivity(:,Inod)+neqx,connectivity(:,Jnod),d2d1(:,Inod,Jnod),neq,neq);
				%K=K+sparse(connectivity(:,Inod)+neqx,connectivity(:,Jnod),d1d2(:,Jnod,Inod),neq,neq);
			end
		end
	end

	
	Kuv=(Kuv+Kuv.')/2 ; % I know that the matrix must be symmetric, but numerically this may not be strickly so
	% Note: for numerical verificatin of distributed parameter gradient it is important to
	% not to use the complex conjugate transpose.
	%whos('K')
	
	% Boundary contribution
	
	if CtrlVar.IncludeBoundaryTerm
		[KBoundary,rhsBoundary]=DirichletBoundaryIntegralDiagnostic(coordinates,connectivity,Boundary,nip,h,u,v,AGlen,n,alpha,rho,rhow,g,CtrlVar);
        
        %save TestSave K KBoundary R rhsBoundary
		Kuv=Kuv+KBoundary ; Ruv=Ruv+rhsBoundary;
	end
end



