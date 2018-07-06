function beta2int=CalcBeta2IntFromNodalValues(u,v,C,m,connectivity,nip,GF)
	
	[Nele,nod]=size(connectivity);
	ndim=2; 
	
	unod=reshape(u(connectivity,1),Nele,nod);
	vnod=reshape(v(connectivity,1),Nele,nod);
	Cnod=reshape(C(connectivity,1),Nele,nod);
	
	[points,weights]=sample('triangle',nip,ndim);
	
	beta2int=zeros(Nele,nip);
	
	for Iint=1:nip
		
		
		fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
		
		uint=unod*fun;
		vint=vnod*fun;
		Cint=Cnod*fun;
		
		beta2int(:,Iint) = calcBeta2in2Dint(uint,vint,Cint,m,GF,Iint);
	end
	
end