function M=MassMatrix2D(coordinates,connectivity,nip) % old version, not used anymore
	
	% old version, not used anymore (MassMatrix2D1dof preferred)
	% calculates the mass matrix, ie : int N_p N_q
	
	
	Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
	ndim=2;  
	
	
	[points,weights]=sample('triangle',nip,ndim);
	d1d1=zeros(Nele,nod,nod); 
	
	
	
	for Iint=1:nip
		fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
		[~,detJ]=derivVector(coordinates,connectivity,nip,Iint);
		detJw=detJ*weights(Iint);
		for Inod=1:nod
			for Jnod=1:nod
				
				d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)...
					+fun(Jnod).*fun(Inod).*detJw;

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
	M=(M+M.')/2 ; % I know that the matrix must be symmetric, but numerically this may not be strickly so

	
end



