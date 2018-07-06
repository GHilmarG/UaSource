function M=StiffnessMatrixSx2D(coordinates,connectivity,nip,CtrlVar)
	
	
	% calculates the xx part of the stiffness mattrix. ie :
	%         p_x N_p \p_x N_q
	
	Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
	ndim=2 ; [points,weights]=sample('triangle',nip,ndim);
	d1d1=zeros(Nele,nod,nod);
	
	beta2=CtrlVar.TotalRegBetaSquare1;
	
	for Iint=1:nip
		[Deriv,detJ]=derivVector(coordinates,connectivity,nip,Iint);
		detJw=detJ*weights(Iint);
		for Inod=1:nod
			for Jnod=1:nod
				d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)...
					+(beta2+Deriv(:,1,Inod).*Deriv(:,1,Jnod)).*detJw;
			end
		end
	end
	
	
	Iind=zeros(nod*nod*Nele,1); Jind=zeros(nod*nod*Nele,1);Xval=zeros(nod*nod*Nele,1); istak=0;
	for Inod=1:nod
		
		for Jnod=1:nod
			
			Iind(istak+1:istak+Nele)=connectivity(:,Inod); Jind(istak+1:istak+Nele)=connectivity(:,Jnod); Xval(istak+1:istak+Nele)=d1d1(:,Inod,Jnod);
			istak=istak+Nele;
			
		end
	end
	
	M=sparse(Iind,Jind,Xval,Nnodes,Nnodes);
	M=(M+M.')/2 ;
end



