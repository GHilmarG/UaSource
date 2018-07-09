function [dIdA]=dIdA_EleSteps(u,v,lx,ly,h,connectivity,coordinates,nip,AGlen,n,CtrlVar)
	
    error('not to be used')
    
	% calculates dIdA with deltaA constant within each element
	%
	
	[Nele,nod]=size(connectivity) ; ndim=2;
	
	
	hnod=reshape(h(connectivity,1),Nele,nod);   % Nele x nod
	unod=reshape(u(connectivity,1),Nele,nod);
	vnod=reshape(v(connectivity,1),Nele,nod);
	lxnod=reshape(lx(connectivity,1),Nele,nod);
	lynod=reshape(ly(connectivity,1),Nele,nod);
	AGlennod=reshape(AGlen(connectivity,1),Nele,nod);
	
	
	[points,weights]=sample('triangle',nip,ndim);
	
	
	dIdA=zeros(Nele,1); EleArea=zeros(Nele,1);
	
	[~,~,~,~,~,~,~,e]=...
		calcStrainRatesEtaInt(u,v,coordinates,connectivity,nip,AGlen,n,CtrlVar);
	
	
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
		dEtadA=-real(hint.*AGlenInt.^(-1/n-1).*e(:,Iint).^((1-n)/n)/(2*n));
		
		EleArea=EleArea+detJw;
		
		dIdA=dIdA-dEtadA.*...
			((4*dudx+2*dvdy).*dlxdx+(dudy+dvdx).*dlxdy+(4*dvdy+2*dudx).*dlydy+(dudy+dvdx).*dlydx).*detJw;

		
	end
	
	% clearly the sensitivity to a constant shift within an element is not independent of the size of the element, other factors equal,
	% a given change in A over a large element can be expected to have a larger impact on J than a change in A over a small element.
	% But in a distributed problem where sensitivities are expressed at nodal points, the desired sensitivities are with respect to a unit area,
	% hence the need to divide with the area of the element.
	
	dIdA=dIdA./EleArea;
end



