function [IRegbeta,IRegb,dIRegdbeta,dIRegdb]=ModelNorm(iModelType,beta,beta_prior,beta_error,bModel,b_prior,b_error,coordinates,connectivity,nip)
                 
    % currently only using iModelType=1
    
    Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
    
    ndim=2; %dof=2; neq=dof*Nnodes;
    
    [points,weights]=sample('triangle',nip,ndim);
    funInt=cell(nip); derInt=cell(nip);
    for Iint=1:nip
        funInt{Iint}=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        derInt{Iint}=shape_der(Iint,ndim,nod,points);  % dof x nod : dNj/dXi=[dN1/dx dN2/dx dN3/dx; dN1/dy dN2/dy dN3/dy]
    end
    
    IReg=0; dIRegdbeta=zeros(Nnodes,1);
    Norm=zeros(Nnodes,1);
    
    if iModelType ==1 % || beta-beta_prior ||
        
        
       
        
        IRegbeta=sum(((beta-beta_prior)/beta_error).^2);
        dIRegdbeta=2*(beta-beta_prior)/beta_error^2;
        
        IRegb=sum(((bModel-b_prior)/b_error).^2);
        dIRegdb=2*(bModel-b_prior)/b_error^2;
        
        IReg=IRegbeta+IRegb; 
        
    elseif iModelType ==2 % || beta-beta_prior ||
        for Iele=1:Nele
            con=connectivity(Iele,:);  % nodes of element
            coo=coordinates(con,:) ; % nod x dof =[x1 y1 ; x2 y2 ; x3 y3]
            
            beta_l=beta(connectivity(Iele,:)) ;
            beta_prior_l=beta_prior(connectivity(Iele,:)) ;
            
            g_l=connectivity(Iele,:);
            
            dIRegdbeta_l=zeros(nod,1);
            Norm_l=zeros(nod,1);
            
            for Iint=1:nip                           % loop over integration points
                
                fun=funInt{Iint} ;
                der=derInt{Iint};
                J=der*coo; % (dof x nod) x (nod x dof) = dof x dof
                detJ=det(J);  % det(dof x dof) matrix
                detJw=detJ*weights(Iint);
                
                betaI=beta_l'*fun;
                beta_priorI=beta_prior_l'*fun;
                
                IReg=IReg+(betaI-beta_priorI)^2*detJw; % integral
                
                dIRegdbeta_l=dIRegdbeta_l+2*(betaI-beta_priorI)*fun*detJw;
                Norm_l=Norm_l+fun*detJw;
            end
            
            for i1=1:length(g_l)
                dIRegdbeta(g_l(i1))=dIRegdbeta(g_l(i1))+dIRegdbeta_l(i1);
                Norm(g_l(i1))=Norm(g_l(i1))+Norm_l(i1);
            end
        end
        
         dIRegdbeta=EleAverageInterpolate(dIRegdbeta,coordinates,connectivity);
        
        
    elseif iModelType ==3
        
        
        %          cooTri=zeros(Nele,2);
        %          dIRegdbetaTri=zeros(Nele,1);
        %
        dIRegdbeta=zeros(Nnodes,1);
        
        for Iele=1:Nele
            con=connectivity(Iele,:);  % nodes of element
            coo=coordinates(con,:) ; % nod x dof =[x1 y1 ; x2 y2 ; x3 y3]
            
            beta_l=beta(connectivity(Iele,:)) ;
            g_l=connectivity(Iele,:);
            dIRegdbeta_l=zeros(nod,1);
            
            for Iint=1:nip                           % loop over integration points
                
                fun=funInt{Iint} ;
                der=derInt{Iint};
                
                J=der*coo;    % (dof x nod) x (nod x dof) = dof x dof
                detJ=det(J);  % det(dof x dof) matrix
                detJw=detJ*weights(Iint);
                deriv=J\der;  % (dof x dof) x (dof x nod) = dof x nod
                
                F=(deriv(1,:)*beta_l)^2+(deriv(2,:)*beta_l)^2;
                
                IReg=IReg+F*detJw;
                
                
                F=2*deriv(1,:)*beta_l*deriv(1,:)'+2*deriv(2,:)*beta_l*deriv(2,:)';
                dIRegdbeta_l=dIRegdbeta_l+F*detJw;
                
                
            end
            
            
            for i1=1:length(g_l)
                dIRegdbeta(g_l(i1))=dIRegdbeta(g_l(i1))+dIRegdbeta_l(i1);
            end
            
            
        end
         dIRegdbeta=EleAverageInterpolate(dIRegdbeta,coordinates,connectivity);
    end
    
    IReg=real(IReg);
    dIRegdbeta=real(dIRegdbeta);
    
    
    % calculate average values for each triangle, then interpolate back onto nodes
    
    
   
    
end

