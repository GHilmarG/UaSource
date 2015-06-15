function dIdA=dIdAq(u,v,lx,ly,h,connectivity,coordinates,nip,AGlen,n,CtrlVar)
	
    
    %
    % old test version, not to be used
    %
    
	Nnodes=max(connectivity(:));
	[Nele,nod]=size(connectivity) ; ndim=2;
	
	
	hnod=reshape(h(connectivity,1),Nele,nod);   % Nele x nod
	unod=reshape(u(connectivity,1),Nele,nod);
	vnod=reshape(v(connectivity,1),Nele,nod);
	lxnod=reshape(lx(connectivity,1),Nele,nod);
	lynod=reshape(ly(connectivity,1),Nele,nod);
	AGlennod=reshape(AGlen(connectivity,1),Nele,nod);
	
	
	[points,weights]=sample('triangle',nip,ndim);
	[~,~,~,~,~,~,~,e]=...
		calcStrainRatesEtaInt(u,v,coordinates,connectivity,nip,AGlen,n,CtrlVar);
	
	T=zeros(Nele,nod);
	
	for Iint=1:nip
		
		
		fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
		[Deriv,detJ]=derivVector(coordinates,connectivity,nip,Iint);
		
		hint=hnod*fun;
		AGlenInt=AGlennod*fun;
		
		dudx=zeros(Nele,1); dvdx=zeros(Nele,1); dudy=zeros(Nele,1); dvdy=zeros(Nele,1);
		dlxdx=zeros(Nele,1); dlydx=zeros(Nele,1); dlxdy=zeros(Nele,1); dlydy=zeros(Nele,1);
		
		for Inod=1:nod
			
			dudx=dudx+Deriv(:,1,Inod).*unod(:,Inod);
			dvdx=dvdx+Deriv(:,1,Inod).*vnod(:,Inod);
			dudy=dudy+Deriv(:,2,Inod).*unod(:,Inod);
			dvdy=dvdy+Deriv(:,2,Inod).*vnod(:,Inod);
			
			dlxdx=dlxdx+Deriv(:,1,Inod).*lxnod(:,Inod);
			dlydx=dlydx+Deriv(:,1,Inod).*lynod(:,Inod);
			dlxdy=dlxdy+Deriv(:,2,Inod).*lxnod(:,Inod);
			dlydy=dlydy+Deriv(:,2,Inod).*lynod(:,Inod);
			
		end
		
		detJw=detJ*weights(Iint);
		dEtadA=-real(hint.*AGlenInt.^(-1/n-1).*e(:,Iint).^((1-n)/n))/(2*n);
		
		for Inod=1:nod
			T(:,Inod)=T(:,Inod)...
				-dEtadA.*((4*dudx+2*dvdy).*dlxdx+(dudy+dvdx).*dlxdy+(4*dvdy+2*dudx).*dlydy+(dudy+dvdx).*dlydx).*fun(Inod).*detJw;
		end
	end
	
	dIdA=zeros(Nnodes,1);

	for Inod=1:nod
		dIdA=dIdA+sparse(connectivity(:,Inod),ones(Nele,1),T(:,Inod),Nnodes,1);
	end

	
end



