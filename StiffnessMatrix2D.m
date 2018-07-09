function M=StiffnessMatrix2D(coordinates,connectivity,nip,CtrlVar)
	
	
	% calculates the stiffness mattrix. ie :
	%         p_x N_p \p_x N_q      p_x N_p \p_y N_q
	%         p_y N_p \p_x N_q      p_y N_p \p_y N_q
	
	Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
	ndim=2; dof=2; neq=dof*Nnodes;
	neqx=Nnodes ;
	[points,weights]=sample('triangle',nip,ndim);
	
	d1d1=zeros(Nele,nod,nod); d2d2=zeros(Nele,nod,nod);  d1d2=zeros(Nele,nod,nod); d2d1=zeros(Nele,nod,nod);
	
	
	for Iint=1:nip
		
		[Deriv,detJ]=derivVector(coordinates,connectivity,nip,Iint);
		detJw=detJ*weights(Iint);
		
		for Inod=1:nod
			for Jnod=1:nod
				
				d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)...
					+Deriv(:,1,Inod).*Deriv(:,1,Jnod).*detJw;
				
				d1d2(:,Inod,Jnod)=d1d2(:,Inod,Jnod)...
					+Deriv(:,1,Inod).*Deriv(:,2,Jnod).*detJw;
				
				d2d1(:,Inod,Jnod)=d2d1(:,Inod,Jnod)...
					+Deriv(:,2,Inod).*Deriv(:,1,Jnod).*detJw;
				
				d2d2(:,Inod,Jnod)=d2d2(:,Inod,Jnod)...
					+Deriv(:,2,Inod).*Deriv(:,2,Jnod).*detJw;
				
			end
		end
	end
	
	
	Iind=zeros(nod*nod*Nele*4,1); Jind=zeros(nod*nod*Nele*4,1);Xval=zeros(nod*nod*Nele*4,1); istak=0;
	for Inod=1:nod
		
		for Jnod=1:nod
			
			Iind(istak+1:istak+Nele)=connectivity(:,Inod); Jind(istak+1:istak+Nele)=connectivity(:,Jnod); Xval(istak+1:istak+Nele)=d1d1(:,Inod,Jnod);
			istak=istak+Nele;
			
			Iind(istak+1:istak+Nele)=connectivity(:,Inod)+neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod)+neqx; Xval(istak+1:istak+Nele)=d2d2(:,Inod,Jnod);
			istak=istak+Nele;
			
			Iind(istak+1:istak+Nele)=connectivity(:,Inod); Jind(istak+1:istak+Nele)=connectivity(:,Jnod)+neqx; Xval(istak+1:istak+Nele)=d1d2(:,Inod,Jnod);
			istak=istak+Nele;
			
			Iind(istak+1:istak+Nele)=connectivity(:,Inod)+neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod); Xval(istak+1:istak+Nele)=d2d1(:,Inod,Jnod);
			istak=istak+Nele;
			
		end
	end
	% nzmax=size(unique([Iind Jind],'rows'),1) ; M=sparse(Iind,Jind,Xval,neq,neq,nzmax); not sure why this
	% does not work
	
	M=sparse(Iind,Jind,Xval,neq,neq);
	
	M=(M+M.')/2 ; % I know that the matrix must be symmetric, but numerically this may not be strickly so
	% Note: for numerical verificatin of distributed parameter gradient it is important to
	% not to use the complex conjugate transpose.
	%whos('M')
	
	
end



