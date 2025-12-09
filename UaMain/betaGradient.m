function [betaGrad]=betaGradient(u,v,Lambda,Mu,C,connectivity,coordinates,nip);
  
    disp(' beta gradient ')
        
    Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
    ndim=2; dof=2; neq=dof*Nnodes;
    
    [points,weights]=sample('triangle',nip,ndim);
    funInt=cell(nip); derInt=cell(nip);
    for Iint=1:nip
        funInt{Iint}=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        derInt{Iint}=shape_der(Iint,ndim,nod,points);  % dof x nod : dNj/dXi=[dN1/dx dN2/dx dN3/dx; dN1/dy dN2/dy dN3/dy]
    end
        
    betaGrad=zeros(Nnodes,1);
        
    for Iele=1:Nele
        con=connectivity(Iele,:);  % nodes of element
        coo=coordinates(con,:) ; % nod x dof =[x1 y1 ; x2 y2 ; x3 y3]
        
        u_l=u(connectivity(Iele,:)); v_l=v(connectivity(Iele,:)) ;
        Lambda_l=Lambda(connectivity(Iele,:)); Mu_l=Mu(connectivity(Iele,:));
        C_l=C(con);
        
        g_l=connectivity(Iele,:);
        betaGradI=zeros(nod,1);
        
        for Iint=1:nip                           % loop over integration points
            fun=funInt{Iint} ; der=derInt{Iint};
            J=der*coo; % (dof x nod) x (nod x dof) = dof x dof
            detJ=det(J);  % det(dof x dof) matrix
            detJw=detJ*weights(Iint);
            CI=C_l'*fun ;
            betaI= sqrt(1/CI) ;
            uI=u_l'*fun ;
            vI=v_l'*fun ;
            LambdaI=Lambda_l'*fun;
            MuI=Mu_l'*fun;
            betaGradI=betaGradI-2*(LambdaI*uI+MuI*vI)*betaI*fun*detJw;
            
            
        end
        
        for i1=1:length(g_l)
            betaGrad(g_l(i1))=betaGrad(g_l(i1))+betaGradI(i1);
        end
    end
    
end

