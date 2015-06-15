function [K,R,T,F]=KRTFtest(s,h,u,v,AGlen,n,C,m,coordinates,connectivity,nip,GF,alpha,rho,rhow,g,NLit)
	
	
	% calculates the tangent matrix and right-hand side in a vectorized form
	
	
	
	
	Dvisk=NLit.NRvisk ; % if gradients with respect to visk not to be included set to 0, othewise 1
	Dbeta=NLit.NRbeta2;
	
    
	[beta2,Dbeta2Duu,Dbeta2Dvv,Dbeta2Duv] = calcBeta2in2D(u,v,C,m,GF) ;
	[etaInt,~,~,exx,eyy,exy,Eint]=calcStrainRatesEtaInt(u,v,coordinates,connectivity,nip,AGlen,n);
	
	
	
	Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
	ndim=2; dof=2; neq=dof*Nnodes;
	neqx=Nnodes ;
	
	hnod=reshape(h(connectivity,1),Nele,nod);   % Nele x nod
	snod=reshape(s(connectivity,1),Nele,nod);
	unod=reshape(u(connectivity,1),Nele,nod);
	vnod=reshape(v(connectivity,1),Nele,nod);
	
	Cnod=reshape(C(connectivity,1),Nele,nod);
	
	beta2nod=reshape(beta2(connectivity,1),Nele,nod);
	%Dbeta2nod=reshape(Dbeta2(connectivity,1),Nele,nod);
	Dbeta2Duunod=reshape(Dbeta2Duu(connectivity,1),Nele,nod);
	Dbeta2Dvvnod=reshape(Dbeta2Dvv(connectivity,1),Nele,nod);
	Dbeta2Duvnod=reshape(Dbeta2Duv(connectivity,1),Nele,nod);
	
	
	rhog=rho*g; ca=cos(alpha); sa=sin(alpha);
	
	
	[points,weights]=sample('triangle',nip,ndim);
	% get local coordinates and weights
	
	
	
	
	
	
	
	d1d1=zeros(Nele,nod,nod); d2d2=zeros(Nele,nod,nod);  d1d2=zeros(Nele,nod,nod); d2d1=zeros(Nele,nod,nod);
	
	Tx=zeros(Nele,nod);  Ty=zeros(Nele,nod); Fx=zeros(Nele,nod);  Fy=zeros(Nele,nod);
	
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
		
		[beta2int,Dbeta2Duuint,Dbeta2Dvvint,Dbeta2Duvint] = calcBeta2in2Dint(uint,vint,Cint,m,GF,Iint);
		%beta2int=beta2nod*fun;       %beta2int=beta2int.*gfint(:,Iint);
		
		%Dbeta2Duuint=Dbeta2Duunod*fun; %Dbeta2Duint=Dbeta2Duint.*gfint(:,Iint);
		%Dbeta2Dvvint=Dbeta2Dvvnod*fun; %Dbeta2Dvint=Dbeta2Dvint.*gfint(:,Iint);
		%Dbeta2Duvint=Dbeta2Duvnod*fun; %Dbeta2Dvint=Dbeta2Dvint.*gfint(:,Iint);
		
		%Dbeta2Duint=(Dbeta2nod*fun); %.*uint;
		%Dbeta2Dvint=(Dbeta2nod*fun); %.*vint;
		
		etaint=etaInt(:,Iint) ;  % I could consider calculating this here
		
		
		% derivatives at this integration point for all elements
		dsdx=zeros(Nele,1); dhdx=zeros(Nele,1); dsdy=zeros(Nele,1); dhdy=zeros(Nele,1);
		%dudx=zeros(Nele,1); dudy=zeros(Nele,1); dvdx=zeros(Nele,1); dvdy=zeros(Nele,1);
		
		
		
		for Inod=1:nod
			
			dsdx=dsdx+Deriv(:,1,Inod).*snod(:,Inod);
			dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
			dsdy=dsdy+Deriv(:,2,Inod).*snod(:,Inod);
			dhdy=dhdy+Deriv(:,2,Inod).*hnod(:,Inod);
			
			
		end
		
		% this takes longer!
% 		tic
% 		testdsdx=sum(squeeze(Deriv(:,1,:)).*snod,2);
% 		testdhdx=sum(squeeze(Deriv(:,1,:)).*hnod,2);
% 		testdsdy=sum(squeeze(Deriv(:,2,:)).*snod,2);
% 		testdhdy=sum(squeeze(Deriv(:,2,:)).*hnod,2);
% 		toc
		
		
		
		detJw=detJ*weights(Iint);
		
		
		for Inod=1:nod
			for Jnod=1:nod
				
				
				d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)...
					+(4*hint.*etaint.*Deriv(:,1,Inod).*Deriv(:,1,Jnod)...
					+hint.*etaint.*Deriv(:,2,Inod).*Deriv(:,2,Jnod)...
					+beta2int.*fun(Jnod).*fun(Inod)...
					+Dbeta.*Dbeta2Duuint.*fun(Jnod).*fun(Inod)).*detJw;  % beta derivative	
%					+Dbeta*uint.*Dbeta2Duint.*fun(Jnod).*fun(Inod)).*detJw;  % beta derivative (slower
%					convergence!)
				
				% was
				%+Dbeta*uint.*Dbeta2Duint*fun(Jnod).*fun(Inod)).*detJw;  % beta derivative
				% changed this on 30 July 2010 (is this correct?)
				
				d2d2(:,Inod,Jnod)=d2d2(:,Inod,Jnod)...
					+(4*hint.*etaint.*Deriv(:,2,Inod).*Deriv(:,2,Jnod)...
					+hint.*etaint.*Deriv(:,1,Inod).*Deriv(:,1,Jnod)...
					+beta2int.*fun(Jnod).*fun(Inod)...
					+Dbeta.*Dbeta2Dvvint.*fun(Jnod).*fun(Inod)).*detJw ;
				%+Dbeta*vint.*Dbeta2Dvint.*fun(Jnod).*fun(Inod)).*detJw ;
				
				
				d1d2(:,Inod,Jnod)=d1d2(:,Inod,Jnod)...
					+(etaint.*hint.*(2*Deriv(:,1,Inod).*Deriv(:,2,Jnod)+Deriv(:,2,Inod).*Deriv(:,1,Jnod))...
					+Dbeta.*Dbeta2Duvint.*fun(Jnod).*fun(Inod)).*detJw;    % beta derivative
				%+Dbeta*uint.*Dbeta2Dvint.*fun(Jnod).*fun(Inod)).*detJw;    % beta derivative
				
				d2d1(:,Inod,Jnod)=d2d1(:,Inod,Jnod)...
					+(etaint.*hint.*(2*Deriv(:,2,Inod).*Deriv(:,1,Jnod)+Deriv(:,1,Inod).*Deriv(:,2,Jnod))...
					+Dbeta.*Dbeta2Duvint*fun(Jnod).*fun(Inod)).*detJw;    % beta derivative
			%	+Dbeta*vint.*Dbeta2Duint*fun(Jnod).*fun(Inod)).*detJw;    % beta derivative
				
				
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
			
		
			
			% is it not better to multiply all nonlinear terms at the nodes
			% and then evaluate that product at the integration points rather than evaluating individual terms
			% at the integration points and then multiply?
			
			t1=-rhog*hint.*((dsdx-(1-rho/rhow).*dhdx)*ca-sa).*fun(Inod);
			t2=0.5*ca*rhog.*(1-rho/rhow).*hint.^2.*Deriv(:,1,Inod);
			t3=-hint.*etaint.*(4*exx(:,Iint)+2*eyy(:,Iint)).*Deriv(:,1,Inod);
			t4=-hint.*etaint.*2.*exy(:,Iint).*Deriv(:,2,Inod);
			t5=-beta2int.*uint.*fun(Inod);
			
			Tx(:,Inod)=Tx(:,Inod)-(t3+t4+t5).*detJw;
			Fx(:,Inod)=Fx(:,Inod)+(t1+t2).*detJw;
			
			
			
			t1=-rhog.*hint.*(dsdy-(1-rho/rhow).*dhdy).*ca.*fun(Inod);
			t2=0.5*ca*rhog.*(1-rho/rhow).*hint.^2.*Deriv(:,2,Inod);
			t3=-hint.*etaint.*(4*eyy(:,Iint)+2*exx(:,Iint)).*Deriv(:,2,Inod);
			t4=-hint.*etaint.*2.*exy(:,Iint).*Deriv(:,1,Inod);
			t5=-beta2int.*vint.*fun(Inod);
			
			Ty(:,Inod)=Ty(:,Inod)-(t3+t4+t5).*detJw;
			Fy(:,Inod)=Fy(:,Inod)+(t1+t2).*detJw;
			
			
			
			
			
		end
	end
	
	
	
	
	% assemble right-hand side
	
	T=sparse(neq,1); F=sparse(neq,1);
	
	for Inod=1:nod
		
		
		T=T+sparse(connectivity(:,Inod),ones(Nele,1),Tx(:,Inod),neq,1);
		T=T+sparse(connectivity(:,Inod)+neqx,ones(Nele,1),Ty(:,Inod),neq,1);
		
		F=F+sparse(connectivity(:,Inod),ones(Nele,1),Fx(:,Inod),neq,1);
		F=F+sparse(connectivity(:,Inod)+neqx,ones(Nele,1),Fy(:,Inod),neq,1);
	end
	
	R=T-F;
	
	K=spalloc(neq,neq,neq*25); % not sure how much space is needed but this seems reasonable
	%K=sparse(neq,neq);
	%whos('K')
	% assemble matrix
	iSparse=1;  %	faster
	if iSparse==1
		% uses the sparse function less often
		
		
		Iind=zeros(nod*Nele*4,1); Jind=zeros(nod*Nele*4,1);Xval=zeros(nod*Nele*4,1);
		for Inod=1:nod
			istak=0;
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
			K=K+sparse(Iind,Jind,Xval,neq,neq);
			
			
		end
		
	else
		% creates the sparse matrix in steps, requires no extra
		for Inod=1:nod
			for Jnod=1:nod
				K=K+sparse(connectivity(:,Inod),connectivity(:,Jnod),d1d1(:,Inod,Jnod),neq,neq);
				K=K+sparse(connectivity(:,Inod)+neqx,connectivity(:,Jnod)+neqx,d2d2(:,Inod,Jnod),neq,neq);
				K=K+sparse(connectivity(:,Inod),connectivity(:,Jnod)+neqx,d1d2(:,Inod,Jnod),neq,neq);
				K=K+sparse(connectivity(:,Inod)+neqx,connectivity(:,Jnod),d2d1(:,Inod,Jnod),neq,neq);
				%K=K+sparse(connectivity(:,Inod)+neqx,connectivity(:,Jnod),d1d2(:,Jnod,Inod),neq,neq);
				
			end
		end
		
	end
	
	
	
	K=(K+K.')/2 ; % I know that the matrix must be symmetric, but numerically this may not be strickly so
	% Note: for numerical verificatin of distributed parameter gradient it is important to
	% not to use the complex conjugate transpose.
	%whos('K')
end



