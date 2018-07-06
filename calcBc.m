function [Bc]=calcBc(S,B,h,u,v,C,m,coordinates,connectivity,nip,rho,rhow,CtrlVar)
		
% calculates the matrix \p r /p C
%
% Bc is 2*Nnodes \times Nnodes
%

	Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
	ndim=2; dof=2; 
	neqx=Nnodes ;
	
	hnod=reshape(h(connectivity,1),Nele,nod);   % Nele x nod
	
	unod=reshape(u(connectivity,1),Nele,nod);
	vnod=reshape(v(connectivity,1),Nele,nod);
	Cnod=reshape(C(connectivity,1),Nele,nod);
	Snod=reshape(S(connectivity,1),Nele,nod);
	Bnod=reshape(B(connectivity,1),Nele,nod);

	
	[points,weights]=sample('triangle',nip,ndim);
	
	
	kxC=zeros(Nele,nod,nod); kyC=zeros(Nele,nod,nod);		
	
	for Iint=1:nip
		
		
		fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
		
		[~,detJ]=derivVector(coordinates,connectivity,nip,Iint);
		% Deriv : Nele x dof x nod
		%  detJ : Nele
		
		% values at integration this point
		hint=hnod*fun;
		%sint=snod*fun;
		uint=unod*fun;
		vint=vnod*fun;
		Cint=Cnod*fun;
		Sint=Snod*fun;
		Bint=Bnod*fun;
		
		hfint=(Sint-Bint)*rhow/rho;
		
		kH=CtrlVar.kH;
		Heint = HeavisideApprox(kH,hint-hfint);
		detJw=detJ*weights(Iint);
		
		
		for Inod=1:nod
			for Jnod=1:nod
				
				Ctemp= (1/m)*Heint.*Cint.^(-1/m-1).*(sqrt(uint.*uint+vint.*vint+CtrlVar.SpeedZero^2)).^(1/m-1) ;
				
				
				kxC(:,Inod,Jnod)=kxC(:,Inod,Jnod)+Ctemp.*uint.*fun(Jnod).*fun(Inod).*detJw;
				kyC(:,Inod,Jnod)=kyC(:,Inod,Jnod)+Ctemp.*vint.*fun(Jnod).*fun(Inod).*detJw;  
				
			end
			
		end
	end
	% assemble matrix

	Bc=sparse(dof*Nnodes,Nnodes);
	Iind=zeros(nod*Nele*2,1); Jind=zeros(nod*Nele*2,1);Xval=zeros(nod*Nele*2,1);
	for Inod=1:nod
		istak=0;
		for Jnod=1:nod
			
			Iind(istak+1:istak+Nele)=connectivity(:,Inod); Jind(istak+1:istak+Nele)=connectivity(:,Jnod); Xval(istak+1:istak+Nele)=kxC(:,Inod,Jnod);
			istak=istak+Nele;
			
			Iind(istak+1:istak+Nele)=connectivity(:,Inod)+neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod); Xval(istak+1:istak+Nele)=kyC(:,Inod,Jnod);
			istak=istak+Nele;
			
		end
		Bc=Bc+sparse(Iind,Jind,Xval,2*Nnodes,Nnodes);
	end

end



