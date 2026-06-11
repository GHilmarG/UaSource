function [K,rh]=DirichletBoundaryIntegralDiagnostic(coordinates,connectivity,Boundary,nip,h,u,v,AGlen,n,alpha,rho,rhow,g,CtrlVar)
	
	%%
	%clc ; clear all 
	%load TestSAve coordinates connectivity Boundary h etaInt alpha rho rhow g CtrlVar nip
	%error('dfas')
	% Need to integrate (h eta (4 dudx + 2 dvdy) n_x + h eta (dvdx + dudy) + varrho g h^2 cos(alpha) n_x/2 ) N_p
	% along all edges where I do not impose the natural boundary condition
	%
	%  Boundary.ElementsBCu{I} is a list of elements for which Dirichlet is defined along edge I=1:3
	%
	% loop over each edge
	%   loop over elements
	%
	%    -calculate position s_i and weights of integration points along a 1d line for ndim=1 and nip=nod+1 (at least)
	%
	%   -determine corresponding points in the 2d (eta,xi) plane, so for example along edge 1
	%    well have (s_i,0), along edge 2 it is (0,s_i), along edge 3 it is (1-s ,s )
	%
	%   -calculate integrand just as done in evaluation the stiffness and the mass matrix
	%
	%   -sum over integration points with corresponding weights, as usuall
	%
	Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
	dof=2; neq=dof*Nnodes; ndim=2; 
	neqx=Nnodes ;neqy=Nnodes;
	
	
	rh=zeros(neq,1);
	
	
	
	icount=0;
	for iEdge=1:numel(Boundary.Edge)
		icount=icount+numel(union(Boundary.ElementsBCu{iEdge},Boundary.ElementsBCv{iEdge}));
	end
	
	N=4*nod*nod*icount; Iind=zeros(N,1) ; Jind=Iind ; Xval=Iind;
	
	varrhog=g*rho.*(1-rho/rhow); ca=cos(alpha);
	
	
	
	
	switch CtrlVar.TriNodes   
		case 3 % 1 exact for linear variatoin
			nipEdge=1; 
		case 6   % 3 exact for second degree polynomials
			nipEdge=4; 
		case 10 % mini
			nipEdge=4; 
		otherwise
			error(' case not recognised, TriNodes value incorrect')
	end
	
	istak=0;
	
	
	for iEdge=1:numel(Boundary.Edge)  % loop over edges
		
		[points,weights]=sampleEdge('line',nipEdge,ndim,iEdge);
	   % get local coordinates and weights for gamma along the 1d line from -1 to 1
	
	
		for Iele=union(Boundary.ElementsBCu{iEdge},Boundary.ElementsBCv{iEdge})
			
			% gather local quantities from global arrays
			% note the nodal numbering is clockwise!
			con=connectivity(Iele,:);  % nodes of edge of the elemetn
			coo=coordinates(con,:) ; % nod x dof =[x1 y1 ; x2 y2 ; x3 y3]
			
			h_l=h(con);  u_l=u(con); v_l=v(con); 
			gx_l=con; gy_l=neqx+con;
			
			c11=zeros(nod,nod) ; c12=c11 ; c21=c11 ;c22=c11;
			b1=zeros(nod,1) ; b2=b1 ;
			
			for Iint=1:nipEdge                           % loop over integration points
				
				fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
				der=shape_der(Iint,ndim,nod,points);  % dof x nod : dNj/dXi=[dN1/dx dN2/dx dN3/dx; dN1/dy dN2/dy dN3/dy]
				
				
				J=der*coo; % (dof x nod) x (nod x dof) = dof x dof
				%detJ=det(J);
				%iJ=inv(J); % dof x dof matrix
				deriv=J\der; % (dof x dof) x (dof x nod) = dof x nod
				
				hint=h_l'*fun ;
				
				
				
				nxdGamma(1)=der(1,:)*coo(:,2); % (1,nod) x (nod,1)=scalar at each integration point
				nxdGamma(2)=-der(2,:)*coo(:,2);
				nxdGamma(3)=-(der(1,:)-der(2,:))*coo(:,2);
                
                nydGamma(1)=-der(1,:)*coo(:,1);
                nydGamma(2)=der(2,:)*coo(:,1);
                nydGamma(3)=(der(1,:)-der(2,:))*coo(:,1);
                
                %fprintf('\n \n Iele %i edge %i Iint %i nx %g ny %g \n',Iele,iEdge,Iint,nxdGamma(iEdge),nydGamma(iEdge))
                %connectivity(Iele,:)
                
                
                
                
                exx=deriv(1,:)*u_l;
                eyy=deriv(2,:)*v_l;
                exy=0.5*(deriv(2,:)*u_l+deriv(1,:)*v_l);
                e=sqrt(CtrlVar.EpsZero+exx^2+eyy^2+exx*eyy+exy^2); 
                etaint=real(0.5*AGlen^(-1/n)*e^((1-n)/n));
                
                
    
				c11=c11+etaint*hint*(4*fun*deriv(1,:)*nxdGamma(iEdge)+fun*deriv(2,:)*nydGamma(iEdge))*weights(Iint);
				c12=c12+etaint*hint*(2*fun*deriv(2,:)*nxdGamma(iEdge)+fun*deriv(1,:)*nydGamma(iEdge))*weights(Iint);
				c21=c21+etaint*hint*(2*fun*deriv(1,:)*nydGamma(iEdge)+fun*deriv(2,:)*nxdGamma(iEdge))*weights(Iint);
				c22=c22+etaint*hint*(4*fun*deriv(2,:)*nydGamma(iEdge)+fun*deriv(1,:)*nxdGamma(iEdge))*weights(Iint);
				
				
				b1=b1-0.5*varrhog*ca*hint.*hint*fun*weights(Iint)*nxdGamma(iEdge);
				b2=b2-0.5*varrhog*ca*hint.*hint*fun*weights(Iint)*nydGamma(iEdge);
				
				
			end % integration points
			
			
			for i1=1:length(gx_l)  ;
				for i2=1:length(gx_l)
					istak=istak+1; Iind(istak)=gx_l(i1); Jind(istak)=gx_l(i2); Xval(istak)=c11(i1,i2);
				end
			end
			
			
			for i1=1:length(gy_l)  ;
				for i2=1:length(gy_l)
					istak=istak+1; Iind(istak)=gy_l(i1); Jind(istak)=gy_l(i2); Xval(istak)=c22(i1,i2);
				end
			end
			
			
			for i1=1:length(gx_l)  ;
				for i2=1:length(gy_l)
					istak=istak+1; Iind(istak)=gx_l(i1); Jind(istak)=gy_l(i2); Xval(istak)=c12(i1,i2);
				end
			end
			
			for i1=1:length(gy_l)  ;
				for i2=1:length(gx_l)
					istak=istak+1; Iind(istak)=gy_l(i1); Jind(istak)=gx_l(i2); Xval(istak)=c21(i1,i2);
				end
			end
			
			
			for i1=1:length(gx_l)
				rh(gx_l(i1))=rh(gx_l(i1))+b1(i1);
				rh(gy_l(i1))=rh(gy_l(i1))+b2(i1);
			end
			
		end  % element loop
	end
	
	
	K=sparse(Iind,Jind,Xval,neqx+neqy,neqx+neqy);
	
	%K=(K+K')/2 ; % I know that the matrix must be symmetric, but numerically this may not be strickly so
	
	%%
	
	
end



