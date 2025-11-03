function [kv,rh]=Next2DAssembleMatrix(dt,h0,u0,v0,a0,u1,v1,a1,coordinates,connectivity,nip,CtrlVar)
	
%
%  dh/dt + p_x (u h ) + \p_y (v h) = a  
%
% Theta method: \Delta h / Delta t = \Theta d h_1/dt + (1-\Theta) d h_0/dt 
%
% 0 : Start of time step
% 1 : End of time step
%
% d h_0/dt= a0 - d (u0 h0)/dx - d (v0 h0) /dy 
%
% This is a linear equation with respect to h
%
% 
	
	Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
	ndim=2; dof=1; neq=dof*Nnodes;
	
	theta=CtrlVar.theta;
	
	h0nod=reshape(h0(connectivity,1),Nele,nod);
	u0nod=reshape(u0(connectivity,1),Nele,nod);   % Nele x nod
	u1nod=reshape(u1(connectivity,1),Nele,nod);
	v0nod=reshape(v0(connectivity,1),Nele,nod);   % Nele x nod
	v1nod=reshape(v1(connectivity,1),Nele,nod);
	a0nod=reshape(a0(connectivity,1),Nele,nod);
	a1nod=reshape(a1(connectivity,1),Nele,nod);
	
	
	[points,weights]=sample('triangle',nip,ndim);
	
	kv=sparse(neq,neq);
	d1d1=zeros(Nele,nod,nod);
	b1=zeros(Nele,nod);
	
	% vector over all elements for each integration point
	for Iint=1:nip
		
		
		
		fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
		[Deriv,detJ]=derivVector(coordinates,connectivity,nip,Iint);
		
		% Deriv : Nele x dof x nod
		%  detJ : Nele
		
		% values at integration point
		
		h0int=h0nod*fun;
		u0int=u0nod*fun;
		v0int=v0nod*fun;
		a0int=a0nod*fun;
		
		u1int=u1nod*fun;
		v1int=v1nod*fun;
		a1int=a1nod*fun;
		
		du1dx=zeros(Nele,1); du0dx=zeros(Nele,1); dh0dx=zeros(Nele,1);
		dv1dy=zeros(Nele,1); dv0dy=zeros(Nele,1); dh0dy=zeros(Nele,1);
		
		% derivatives at one integration point for all elements 
		for Inod=1:nod
			du1dx=du1dx+Deriv(:,1,Inod).*u1nod(:,Inod);
			du0dx=du0dx+Deriv(:,1,Inod).*u0nod(:,Inod);
			
			dv1dy=dv1dy+Deriv(:,2,Inod).*v1nod(:,Inod);
			dv0dy=dv0dy+Deriv(:,2,Inod).*v0nod(:,Inod);
			
			dh0dx=dh0dx+Deriv(:,1,Inod).*h0nod(:,Inod);
			dh0dy=dh0dy+Deriv(:,2,Inod).*h0nod(:,Inod);
			
		end
		
		detJw=detJ*weights(Iint);
		
		% dt theta ( d(u1 h1)/dx    + d(v1 h1)/dy) + h1=
		%  h0+dt { (1-theta) a0+theta a1-(1-theta) (d(u0 h0)/dx+d(v0 h0)/dy}
		
		
		for Inod=1:nod
			for Jnod=1:nod
				
				h1term=fun(Jnod).*fun(Inod).*detJw;
				
				hdxu1=dt*theta*du1dx.*fun(Jnod).*fun(Inod).*detJw;
				udxh1=dt*theta*u1int.*Deriv(:,1,Jnod).*fun(Inod).*detJw;
				
				hdyv1=dt*theta*dv1dy.*fun(Jnod).*fun(Inod).*detJw;
				vdyh1=dt*theta*v1int.*Deriv(:,2,Jnod).*fun(Inod).*detJw;
				
				d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+h1term+hdxu1+udxh1+hdyv1+vdyh1;
				
			end
			
			h0term=h0int.*fun(Inod).*detJw;  % this is the h term, ie not \Delta h term (because the system is linear in h no need to write it in incremental form)
			a0term=dt*(1-theta)*a0int.*fun(Inod).*detJw;
			a1term=dt*theta*a1int.*fun(Inod).*detJw;
			
			hdxu0=-dt*(1-theta)*du0dx.*h0int.*fun(Inod).*detJw;
			udxh0=-dt*(1-theta)*dh0dx.*u0int.*fun(Inod).*detJw;
			
			hdyv0=-dt*(1-theta)*dv0dy.*h0int.*fun(Inod).*detJw;
			vdyh0=-dt*(1-theta)*dh0dy.*v0int.*fun(Inod).*detJw;
						
			b1(:,Inod)=b1(:,Inod)+h0term+a0term+a1term+hdxu0+udxh0+hdyv0+vdyh0;
			
		end
	end
	
	% assemble right-hand side
	
	rh=sparseUA(neq,1);
	for Inod=1:nod
		rh=rh+sparseUA(connectivity(:,Inod),ones(Nele,1),b1(:,Inod),neq,1);
	end
	
	
	for Inod=1:nod
		for Jnod=1:nod
			kv=kv+sparseUA(connectivity(:,Inod),connectivity(:,Jnod),d1d1(:,Inod,Jnod),neq,neq);
		end
	end
	
	
end