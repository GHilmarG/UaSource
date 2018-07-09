function Prod=FE_f_times_l(f,l,coordinates,connectivity,nip)
    % dIdC=Bc'*l;
    % calculates right-hand side in a vectorized form

    Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
    ndim=2;
	
    
	fx=f(1:Nnodes) ; fy=f(Nnodes+1:2*Nnodes);
	lx=l(1:Nnodes) ; ly=l(Nnodes+1:2*Nnodes);
	
    fxnod=reshape(fx(connectivity,1),Nele,nod);   % Nele x nod
    fynod=reshape(fy(connectivity,1),Nele,nod);
	
	lxnod=reshape(lx(connectivity,1),Nele,nod);
	lynod=reshape(ly(connectivity,1),Nele,nod);

    [points,weights]=sample('triangle',nip,ndim);
    Prod=zeros(Nele,1);
    
    for Iint=1:nip
        
        
        fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        [~,detJ]=derivVector(coordinates,connectivity,nip,Iint);
        
        fxint=fxnod*fun;  fyint=fynod*fun;
		lxint=lxnod*fun;  lyint=lynod*fun;
        detJw=detJ*weights(Iint);
        
        
        for Inod=1:nod
            Prod=Prod+(fxint.*lxint+fyint.*lyint).*detJw;
        end
    end
    
    Prod=sum(Prod);
end



