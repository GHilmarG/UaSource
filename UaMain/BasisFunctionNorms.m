function Nnorm=BasisFunctionNorms(connectivity,coordinates,nip)
    
    disp(' Calculating norm of basis functions ')
    
    Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
    ndim=2; dof=2; neq=dof*Nnodes;
    
    [points,weights]=sample('triangle',nip,ndim);
    funInt=cell(nip); derInt=cell(nip);
   
    for Iint=1:nip
        shape_fun(Iint,ndim,nod,points);
        shape_der(Iint,ndim,nod,points);
        funInt{Iint}=shape_fun(Iint,ndim,nod,points) ;
        derInt{Iint}=shape_der(Iint,ndim,nod,points);
    end
    
    NkNk=zeros(Nnodes,1);
    
    for Iele=1:Nele
        con=connectivity(Iele,:);  % nodes of element
        coo=coordinates(con,:) ; % nod x dof =[x1 y1 ; x2 y2 ; x3 y3]
        g_l=connectivity(Iele,:);
        NkNkEle=zeros(nod,1);
        
        for Iint=1:nip                           % loop over integration points
                   
            fun=funInt{Iint} ; der=derInt{Iint};
            J=der*coo; % (dof x nod) x (nod x dof) = dof x dof
            detJ=det(J);  % det(dof x dof) matrix
            detJw=detJ*weights(Iint);
            
            for Inod=1:nod
                NkNkEle(Inod)=NkNkEle(Inod)+fun(Inod).*fun(Inod).*detJw; % jakobian is always negative so I have to change the sign here
				%NkNkEle(Inod)=NkNkEle(Inod)+fun(Inod).*detJw; % jakobian is always negative so I have to change the sign here
            end
            
            
        end
        
        for i1=1:nod
            NkNk(g_l(i1))=NkNk(g_l(i1))+NkNkEle(i1);
        end
    end
   %Nnorm=NkNk;
   Nnorm=sqrt(NkNk);
    
end

