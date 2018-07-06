function [K,R]=KRTFloop(s,S,B,h,u,v,AGlen,n,C,m,coordinates,connectivity,Boundary,nip,alpha,rho,rhow,g,CtrlVar)
	
	
	% calculates the tangent matrix (K) and right-hand side (-R) in a non-vectorized form
	
	% just done for testing purposes, only implemented n=m=1 (although changing this is easy)
	
	%	R=T-F;
	%   F are the `external forces', i.e. right-hand side of the original system
	%   T are the `internal forces'. The equation is considered solved once internal and external forces are
	%   equal to within a given tolerance
	
	
	
	
	if any(h<0) ; error('KRTF: h negative ') ; end
	if any(C<0) ; error('KRTF: C negative ') ; end
	%if ~isreal(C) ; error('KRTF: C not real ') ; end
	
	
	Dvisk=CtrlVar.NRviscosity ; % if gradients with respect to visk not to be included set to 0, othewise 1
	Dbeta=CtrlVar.NRbeta2;
	
	
	Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
	ndim=2; dof=2; neq=dof*Nnodes;
	neqx=Nnodes ; neqy=Nnodes;
	
	N=4*nod*nod*Nele; Iind=zeros(N,1) ; Jind=Iind ; Xval=Iind;
	
	R=zeros(neq,1) ;
	rhog=rho*g; ca=cos(alpha); sa=sin(alpha);
	
	[etaInt,~,~,exx,eyy,exy,Eint]=calcStrainRatesEtaInt(u,v,coordinates,connectivity,nip,AGlen,n,CtrlVar);
	
	[points,weights]=sample('triangle',nip,ndim);
	
	
	funInt=cell(nip); derInt=cell(nip);
	for Iint=1:nip
		funInt{Iint}=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
		derInt{Iint}=shape_der(Iint,ndim,nod,points);  % dof x nod : dNj/dXi=[dN1/dx dN2/dx dN3/dx; dN1/dy dN2/dy dN3/dy]
	end
	
	
	% [ needed for testing purposes
	hnod=reshape(h(connectivity,1),Nele,nod);   % Nele x nod
	snod=reshape(s(connectivity,1),Nele,nod);
	unod=reshape(u(connectivity,1),Nele,nod);
	vnod=reshape(v(connectivity,1),Nele,nod);
	Cnod=reshape(C(connectivity,1),Nele,nod);
	Snod=reshape(S(connectivity,1),Nele,nod);
	Bnod=reshape(B(connectivity,1),Nele,nod);
	% ]
	
	
	istak=0;
	
	for Iele=1:Nele
		
		con=connectivity(Iele,:);  % nodes of element
		coo=coordinates(con,:) ; % nod x dof =[x1 y1 ; x2 y2 ; x3 y3]
		h_l=h(con); s_l=s(con); u_l=u(con); v_l=v(con) ; C_l=C(con); B_l=B(con) ; S_l=S(con);
		gx_l=con; gy_l=neqx+con;
		
		d1d1=zeros(nod,nod); d2d2=zeros(nod,nod);  d1d2=zeros(nod,nod); d2d1=zeros(nod,nod);
		Tx=zeros(nod,1);  Ty=zeros(nod,1); Fx=zeros(nod,1);  Fy=zeros(nod,1);
		
		
		for Iint=1:nip
			
			fun=funInt{Iint} ; der=derInt{Iint};
			J=der*coo; % (dof x nod) x (nod x dof) = dof x dof
			detJ=det(J);  % det(dof x dof) matrix
			deriv=J\der; % (dof x dof) x (dof x nod) = dof x nod
			
			[Deriv,detJTest]=derivVector(coordinates,connectivity,nip,Iint);
% 			
% 			test1=norm(reshape(Deriv(Iele,:,:),dof,nod)-deriv)/norm(deriv);
% 			test2=norm(detJ-detJTest(Iele))/norm(detJ);
% 			
% 			if test1> 100*eps
% 				test1
% 				reshape(Deriv(Iele,:,:),dof,nod)-deriv
% 				error(' test1 failed ')
% 			end
% 			if test2> 100*eps ;
% 				detJ
% 				detJTest(Iele)
% 				error(' test2 failed ')
% 			end
			
			
			
			etaint=etaInt(Iele,Iint) ;  % scalar
			hint=h_l'*fun ;
			Cint=C_l'*fun ;
			uint=u_l'*fun ;
			vint=v_l'*fun ;
			Bint=B_l'*fun;
			Sint=S_l'*fun;
			
			hfint=rhow*(Sint-Bint)/rho;
			kH=CtrlVar.kH;
			Heint = HeavisideApprox(kH,hint-hfint);
			
			dsdx=deriv(1,:)*s_l;
			dhdx=deriv(1,:)*h_l;
			dsdy=deriv(2,:)*s_l;
			dhdy=deriv(2,:)*h_l;
			
			dsdxT=0 ; dsdyT=0 ; dhdxT=0; dhdyT=0 ;
			for Inod=1:nod
				
				dsdxT=dsdxT+Deriv(Iele,1,Inod).*snod(Iele,Inod);
				dhdxT=dhdxT+Deriv(Iele,1,Inod).*hnod(Iele,Inod);
				dsdyT=dsdyT+Deriv(Iele,2,Inod).*snod(Iele,Inod);
				dhdyT=dhdyT+Deriv(Iele,2,Inod).*hnod(Iele,Inod);
			end
			
			if abs(dsdx-dsdxT) > 100* eps ; error('dsdx') ; end
			if abs(dsdy-dsdyT) > 100* eps ; error('dsdy') ;end
			if abs(dhdx-dhdxT) > 100* eps ; error('dhdx');end
			if abs(dhdy-dhdyT) > 100* eps ; error('dhdy');end
			
			
			
			
			[beta2int,Dbeta2Duuint,Dbeta2Dvvint,Dbeta2Duvint] = calcBeta2in2Dint(uint,vint,Cint,m,Heint,CtrlVar);
			
			
			detJw=detJ*weights(Iint);
			
			d1d1=d1d1...
				+(4*hint*etaint*deriv(1,:)'*deriv(1,:)...
				+hint*etaint*deriv(2,:)'*deriv(2,:)...
				+beta2int*(fun*fun')...
				+Dbeta*Dbeta2Duuint*(fun*fun'))*detJw;  % beta derivative
			
			
			d2d2=d2d2...
				+(4*hint*etaint*deriv(2,:)'*deriv(2,:)...
				+hint*etaint*deriv(1,:)'*deriv(1,:)...
				+beta2int*(fun*fun')...
				+Dbeta*Dbeta2Dvvint*(fun*fun'))*detJw ;
			
			d1d2=d1d2...
				+(etaint*hint*(2*deriv(1,:)'*deriv(2,:)+deriv(2,:)'*deriv(1,:))...
				+Dbeta*Dbeta2Duvint*(fun*fun'))*detJw;    % beta derivative
			
			
			d2d1=d2d1...
				+(etaint*hint*(2*deriv(2,:)'*deriv(1,:)+deriv(1,:)'*deriv(2,:))...
				+Dbeta*Dbeta2Duvint*(fun*fun'))*detJw;    % beta derivative
			
			
			Deu=Eint(Iele,Iint)*((2*exx(Iele,Iint)+eyy(Iele,Iint))*deriv(1,:)'+exy(Iele,Iint)*deriv(2,:)');
			Dev=Eint(Iele,Iint)*((2*eyy(Iele,Iint)+exx(Iele,Iint))*deriv(2,:)'+exy(Iele,Iint)*deriv(1,:)');
			
			% E11=h Deu (4 p_x u + 2 p_y v)   + h Deu  ( p_x v + p_y u) p_y N_p
			
			E11=  hint*(4*exx(Iele,Iint)+2*eyy(Iele,Iint))*Deu*deriv(1,:)...
				+2*hint*exy(Iele,Iint)*Deu*deriv(2,:);
			
			
			E12=  hint*(4*exx(Iele,Iint)+2*eyy(Iele,Iint))*Dev*deriv(1,:)...
				+2*hint*exy(Iele,Iint)*Dev*deriv(2,:);
			
			
			
			E22=  hint*(4*eyy(Iele,Iint)+2*exx(Iele,Iint))*Dev*deriv(2,:)...
				+2*hint*exy(Iele,Iint)*Dev*deriv(1,:);
			
			
			E21= hint*(4*eyy(Iele,Iint)+2*exx(Iele,Iint))*Deu*deriv(2,:)...
				+2*hint*exy(Iele,Iint)*Deu*deriv(1,:);
			
			
			d1d1=d1d1+Dvisk*E11.*detJw;
			d2d2=d2d2+Dvisk*E22.*detJw;
			d1d2=d1d2+Dvisk*E12.*detJw;
			d2d1=d2d1+Dvisk*E21.*detJw;
			
			
			t1=-rhog*hint*((dsdx-(1-rho/rhow)*dhdx)*ca-sa)*fun;
			t2=0.5*ca*rhog*(1-rho/rhow)*hint.^2*deriv(1,:)';
			t3=-hint*etaint*(4*exx(Iele,Iint)+2*eyy(Iele,Iint))*deriv(1,:)';
			t4=-hint*etaint*2*exy(Iele,Iint)*deriv(2,:)';
			t5=-beta2int*uint*fun;
			
			Tx=Tx-(t3+t4+t5)*detJw;
			Fx=Fx+(t1+t2)*detJw;
			
			
			
			t1=-rhog*hint*(dsdy-(1-rho/rhow)*dhdy)*ca*fun;
			t2=0.5*ca*rhog*(1-rho/rhow)*hint.^2*deriv(2,:)';
			t3=-hint*etaint*(4*eyy(Iele,Iint)+2*exx(Iele,Iint))*deriv(2,:)';
			t4=-hint*etaint*2*exy(Iele,Iint)*deriv(1,:)';
			t5=-beta2int*vint*fun;
			
			Ty=Ty-(t3+t4+t5)*detJw;
			Fy=Fy+(t1+t2)*detJw;
			
			
		end
		
		for i1=1:length(gx_l)  ;
			for i2=1:length(gx_l)
				istak=istak+1; Iind(istak)=gx_l(i1); Jind(istak)=gx_l(i2); Xval(istak)=d1d1(i1,i2);
			end
		end
		
		
		for i1=1:length(gy_l)  ;
			for i2=1:length(gy_l)
				istak=istak+1; Iind(istak)=gy_l(i1); Jind(istak)=gy_l(i2); Xval(istak)=d2d2(i1,i2);
			end
		end
		
		
		for i1=1:length(gx_l)  ;
			for i2=1:length(gy_l)
				istak=istak+1; Iind(istak)=gx_l(i1); Jind(istak)=gy_l(i2); Xval(istak)=d1d2(i1,i2);
			end
		end
		
		for i1=1:length(gy_l)  ;
			for i2=1:length(gx_l)
				istak=istak+1; Iind(istak)=gy_l(i1); Jind(istak)=gx_l(i2); Xval(istak)=d2d1(i1,i2);
			end
		end
		
		
		for i1=1:length(gx_l)
			R(gx_l(i1))=R(gx_l(i1))+Tx(i1)-Fx(i1);
			R(gy_l(i1))=R(gy_l(i1))+Ty(i1)-Fy(i1);
		end
		
	end
	
	K=sparse(Iind,Jind,Xval,neqx+neqy,neqx+neqy);
	
	K=(K+K')/2 ; % I know that the matrix must be symmetric, but numerically this may not be strickly so
	
end


