function [R,T,F]=RTF(s,S,B,h,u,v,AGlen,n,C,m,coordinates,connectivity,Boundary,nip,alpha,rho,rhow,g,CtrlVar)
	

	% calculates right-hand side in a vectorized form

	
	[etaInt,~,~,exx,eyy,exy]=calcStrainRatesEtaInt(u,v,coordinates,connectivity,nip,AGlen,n,CtrlVar);
    
    Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
    ndim=2; dof=2; neq=dof*Nnodes;
    neqx=Nnodes ;
    H=S-B;
    hf=rhow*H./rho ;
    He = HeavisideApprox(CtrlVar.kH,h-hf);  % 1 if grounded, 0 if afloat
    
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
        Bint=Bnod*fun;
        Sint=Snod*fun;
        rhoint=rhonod*fun;
        
        hfint=rhow*(Sint-Bint)./rhoint;
        %Heint=hint > hfint;  % not a good idea
        kH=CtrlVar.kH;
        Heint = HeavisideApprox(kH,hint-hfint);
        [beta2int] = calcBeta2in2Dint(uint,vint,Cint,m,Heint,CtrlVar);
        etaint=etaInt(:,Iint) ;  % I could consider calculating this here
        
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

	% assemble right-hand side
	
	T=spalloc(neq,1,neq); F=spalloc(neq,1,neq);
	
	for Inod=1:nod
		
		
		T=T+sparse(connectivity(:,Inod),ones(Nele,1),Tx(:,Inod),neq,1);
		T=T+sparse(connectivity(:,Inod)+neqx,ones(Nele,1),Ty(:,Inod),neq,1);
		
		F=F+sparse(connectivity(:,Inod),ones(Nele,1),Fx(:,Inod),neq,1);
		F=F+sparse(connectivity(:,Inod)+neqx,ones(Nele,1),Fy(:,Inod),neq,1);
	end
	
	R=T-F;
    
	if CtrlVar.IncludeBoundaryTerm
		
        [~,rhsBoundary]=DirichletBoundaryIntegralDiagnostic(coordinates,connectivity,Boundary,nip,h,u,v,AGlen,n,alpha,rho,rhow,g,CtrlVar);
		 R=R+rhsBoundary;
	end
	
end



