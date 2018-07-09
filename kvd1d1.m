function [kv]=...
        kvd1d1(s,h,u,v,coordinates,connectivity,Nnodes,Nele,nip,nod,etaInt,gfint,AGlen,C,...
        L,Lb,lambdau,n,m,alpha,rho,rhow,g,Itime)
  
  
    %if any(h<0) ; error(' thickness negative ') ; end
    
    beta2=zeros(Nele,nip);
    
    calcStrainRates=1;
    
    ndim=2; dof=2; neq=dof*Nnodes;
    neqx=Nnodes ; neqy=Nnodes;
    
    N=4*nod*nod*Nele; Iind=zeros(N,1) ; Jind=Iind ; Xval=Iind;
    
    % solves the non-linear 2d shallow ice stream equation
    % and shallow shelf equation
    % using the method of finite elements
    %
    
    rhog=rho*g; ca=cos(alpha); sa=sin(alpha);
    taub=rhog*sin(alpha)*mean(h);
    %disp([' taub = ',num2str(taub),' u = ',num2str(mean(C)*taub^m)])
    
    [points,weights]=sample('triangle',nip,ndim);
    % get local coordinates and weights
    
    diff=1e10;
    
    %[--  constants related to non-linear loop. Values used might affect solution
    beta2Imax=1e5; tol=1e-5 ;  etamax=1e15 ; e00=1e-10;
    %]-
    
    
    
    if n==1 && m==1 ;
        iteration_max=1 ; iteration_min=1 ; e00=0;
    else
        disp(' non-linear ')
        iteration_max=50; iteration_min=2 ;
    end
    
    iteration=0;
    
    funInt=cell(nip); derInt=cell(nip);
    for Iint=1:nip
        funInt{Iint}=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        derInt{Iint}=shape_der(Iint,ndim,nod,points);  % dof x nod : dNj/dXi=[dN1/dx dN2/dx dN3/dx; dN1/dy dN2/dy dN3/dy]
    end
    
    
    rh=zeros(neq,1) ;
    
    istak=0;
    % element loop
  
    for Iele=1:Nele
        % gather local quantities from global arrays
        % note the nodal numbering is clockwise!
        con=connectivity(Iele,:);  % nodes of element
        coo=coordinates(con,:) ; % nod x dof =[x1 y1 ; x2 y2 ; x3 y3]
        h_l=h(con); s_l=s(con); u_l=u(con); v_l=v(con) ; C_l=C(con);
        gx_l=con; gy_l=neqx+con;
        
        c11=zeros(nod,nod) ; c12=c11 ; c21=c11 ;c22=c11;
        b1=zeros(nod,1) ; b2=b1 ;
        
        for Iint=1:nip                           % loop over integration points
            
            
            % only dpendent on integration points
            
            %fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
            %der=shape_der(Iint,ndim,nod,points);  % dof x nod : dNj/dXi=[dN1/dx dN2/dx dN3/dx; dN1/dy dN2/dy dN3/dy]
            
            fun=funInt{Iint} ; der=derInt{Iint};
            
            J=der*coo; % (dof x nod) x (nod x dof) = dof x dof
            detJ=det(J);  % det(dof x dof) matrix
            deriv=J\der; % (dof x dof) x (dof x nod) = dof x nod
            etaI=etaInt(Iele,Iint) ;  % scalar
            hI=h_l'*fun ;
            CI=C_l'*fun ;
            uI=u_l'*fun ;
            vI=v_l'*fun ;
            
            
            
            if m~=1
                beta2I= CI^(-1/m)* (sqrt(uI.*uI+vI.*vI))^(1/m-1) ; % scalar (think about taking out of loop, and getting rid of uI and vI
                beta2I(beta2I>beta2Imax)=beta2Imax ;
                
            else
                beta2I= 1/CI ;
            end
            beta2I=beta2I*gfint(Iele,Iint);
            beta2(Iele,Iint)=beta2I;
            
                        
            detJw=detJ*weights(Iint);
            d1d1=deriv(1,:)'*deriv(1,:);
            d2d2=deriv(2,:)'*deriv(2,:);
            
            c11=c11-(etaI*hI*(4*d1d1+d2d2)+ beta2I*(fun*fun'))*detJw;
            c12=c12-(etaI*hI*(2*deriv(1,:)'*deriv(2,:)+deriv(2,:)'*deriv(1,:)))*detJw;
            c21=c21-(etaI*hI*(2*deriv(2,:)'*deriv(1,:)+deriv(1,:)'*deriv(2,:)))*detJw;
            %c21=c12;
             c22=c22-(etaI*hI*(4*d2d2+d1d1)+beta2I*(fun*fun'))*detJw;
            
            %b1=b1+rhog*hI*((deriv(1,:)*s_l)*ca-sa)*fun*detJw;
            t1=rhog*hI*((deriv(1,:)*s_l-(1-rho/rhow)*deriv(1,:)*h_l)*ca-sa)*fun;
            t2=-0.5*rhog*(1-rho/rhow)*hI^2*deriv(1,:)';
            b1=b1+(t1+t2)*detJw;
            
            %b2=b2+rhog*hI*(deriv(2,:)*s_l)*ca*fun*detJw;
            t1=rhog*hI*(deriv(2,:)*s_l-(1-rho/rhow)*deriv(2,:)*h_l)*ca*fun;
            t2=-0.5*rhog*(1-rho/rhow)*hI^2*deriv(2,:)';
            b2=b2+(t1+t2)*detJw;
            
        end % integration points
        
        
        % assemble matrix
        %  (if repeating don't calculate Iind adn Jind
        
        
        
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
    
    
    
    kv=sparse(Iind,Jind,Xval,neqx+neqy,neqx+neqy);
    
    kv=(kv+kv')/2 ; % I know that the matrix must be symmetric, but numerically this may not be strickly so
    
    
    
end



