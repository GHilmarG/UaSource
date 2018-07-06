function M=StiffnessMatrixBlockDiagonal2D(coordinates,connectivity,nip,CtrlVar)
	
	
	% calculates the block diagonal stiffness mattrix. ie :
	%         p_x N_p \p_x N_q              0
	%                 0             p_y N_p \p_y N_q
	
	Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
	dof=2; neq=dof*Nnodes;
	neqx=Nnodes ;	ndim=2 ; [points,weights]=sample('triangle',nip,ndim);
	d1d1=zeros(Nele,nod,nod); d2d2=zeros(Nele,nod,nod); 
	
	beta2=CtrlVar.TotalRegBetaSquare1;
	
	for Iint=1:nip
		
		[Deriv,detJ]=derivVector(coordinates,connectivity,nip,Iint);
		detJw=detJ*weights(Iint);
		
		for Inod=1:nod
			for Jnod=1:nod
				
				d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)...
					+(beta2+Deriv(:,1,Inod).*Deriv(:,1,Jnod)).*detJw;
				
				d2d2(:,Inod,Jnod)=d2d2(:,Inod,Jnod)...
					+(beta2+Deriv(:,2,Inod).*Deriv(:,2,Jnod)).*detJw;
				
			end
		end
	end
	
	
	Iind=zeros(nod*nod*Nele*2,1); Jind=zeros(nod*nod*Nele*2,1);Xval=zeros(nod*nod*Nele*2,1); istak=0;
	for Inod=1:nod
		
		for Jnod=1:nod
			
			Iind(istak+1:istak+Nele)=connectivity(:,Inod); Jind(istak+1:istak+Nele)=connectivity(:,Jnod); Xval(istak+1:istak+Nele)=d1d1(:,Inod,Jnod);
			istak=istak+Nele;
			
			Iind(istak+1:istak+Nele)=connectivity(:,Inod)+neqx; Jind(istak+1:istak+Nele)=connectivity(:,Jnod)+neqx; Xval(istak+1:istak+Nele)=d2d2(:,Inod,Jnod);
			istak=istak+Nele;
			
		end
	end
	% nzmax=size(unique([Iind Jind],'rows'),1) ; M=sparse(Iind,Jind,Xval,neq,neq,nzmax); not sure why this
	% does not work
	
	M=sparse(Iind,Jind,Xval,neq,neq);
	
	M=(M+M.')/2 ; 
end



