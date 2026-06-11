function [bGrad]=bGradient(u,v,s,Lambda,Mu,etaInt,alpha,rho,g,connectivity,coordinates,nip);
    
    % an atempt at calculatin the gradiet of J with respect to b, almost certainly not correct !!!
    save bGrad u v s Lambda Mu etaInt alpha rho g connectivity coordinates nip 
    
    
   
    
    Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
    ndim=2; dof=2; neq=dof*Nnodes; 
    
    [points,weights]=sample('triangle',nip,ndim);
    funInt=cell(nip); derInt=cell(nip);
    for Iint=1:nip
        funInt{Iint}=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        derInt{Iint}=shape_der(Iint,ndim,nod,points);  % dof x nod : dNj/dXi=[dN1/dx dN2/dx dN3/dx; dN1/dy dN2/dy dN3/dy]
    end
    
    bGrad=zeros(Nnodes,1);
     
    
    for Iele=1:Nele
        con=connectivity(Iele,:);  % nodes of element
        coo=coordinates(con,:) ; % nod x dof =[x1 y1 ; x2 y2 ; x3 y3]
        
        u_l=u(connectivity(Iele,:)); v_l=v(connectivity(Iele,:)) ;
        Lambda_l=Lambda(connectivity(Iele,:)); Mu_l=Mu(connectivity(Iele,:));
        
        s_l=s(connectivity(Iele,:));
        
        g_l=connectivity(Iele,:);
        bGradI=zeros(nod,1);
         
        
        for Iint=1:nip                           % loop over integration points
            fun=funInt{Iint} ; der=derInt{Iint};
            J=der*coo; % (dof x nod) x (nod x dof) = dof x dof
            detJ=det(J);  % det(dof x dof) matrix
            detJw=detJ*weights(Iint);
            deriv=J\der; % (dof x dof) x (dof x nod) = dof x nod
            
            dudxI=deriv(1,:)*u_l;
            dudyI=deriv(2,:)*u_l;
            dvdxI=deriv(1,:)*v_l;
            dvdyI=deriv(2,:)*v_l;
            dLdxI=deriv(1,:)*Lambda_l;
            dLdyI=deriv(2,:)*Lambda_l;
            dMdxI=deriv(1,:)*Mu_l;
            dMdyI=deriv(2,:)*Mu_l;
            etaI=etaInt(Iele,Iint);
            dsdxI=s_l*deriv(1,:);
            dsdyI=s_l*deriv(2,:);
            
            t1=etaI*((4*dudxI+2*dvdyI)*dLdxI+(dvdxI+dudyI)*dLdyI)+rho*g*(dsdxI*cos(alpha)-sin(alpha));
            t2=etaI*((4*dvdyI+2*dudxI)*dMdyI+(dudyI+dvdxI)*dMdxI)+rho*g*dsdyI*cos(alpha);
            
    
            
           % t1=rho*g*(dsdxI*cos(alpha)-sin(alpha));
           % t2=rho*g*dsdyI*cos(alpha);
           
            
            bGradI=bGradI+(t1+t2)*fun*detJw;
            
        end
        
        for i1=1:length(g_l)
            bGrad(g_l(i1))=bGrad(g_l(i1))+bGradI(i1);
         
        end
    end
    
   
    
    
    
    