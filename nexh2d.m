function [h1]=...
        nexh2d(dt,h0,u0,v0,a0,u1,v1,a1,coordinates,connectivity,Nnodes,Nele,nip,nod,L,Lrhs,theta)
    
    
    
    
    ndim=2; neq=Nnodes ;
    
    
    [points,weights]=sample('triangle',nip,ndim);
    % get local coordinates and weights
    
    rh=zeros(neq,1) ; kv=zeros(neq,neq)  ;
    
    
    % element loop
    for Iele=1:Nele
        % gather local quantities from global arrays
        % note the nodal numbering is clockwise!
        con=connectivity(Iele,:);  % nodes of element
        coo=coordinates(con,:) ; % nod x dof =[x1 y1 ; x2 y2 ; x3 y3]
        h0_l=h0(con); u0_l=u0(con); v0_l=v0(con) ; a0_l=a0(con) ;
        u1_l=u1(con); v1_l=v1(con) ; a1_l=a1(con) ;
        
        km=zeros(nod,nod); rhs=zeros(nod,1);
        
        g_l=con;
        
        for Iint=1:nip                           % loop over integration points
            
            fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
            der=shape_der(Iint,ndim,nod,points);  % dof x nod : dNj/dXi=[dN1/dx dN2/dx dN3/dx; dN1/dy dN2/dy dN3/dy]
            
            J=der*coo; % (dof x nod) x (nod x dof) = dof x dof
            detJ=det(J); iJ=inv(J); % dof x dof matrix
            deriv=iJ*der; % (dof x dof) x (dof x nod) = dof x nod
            
            h0I=h0_l'*fun ; u0I=u0_l'*fun ; v0I=v0_l'*fun ; a0I=a0_l'*fun;
            u1I=u1_l'*fun ; v1I=v1_l'*fun ; a1I=a1_l'*fun;
            
            % h1+dt theta (d(u1 h1)/dx + d(v1 h1)/dy)=h0+dt theta a1 + dt (1-theta) a0
            %    -dt (1-theta) (d(u0 h0)/dx + d(v0 h0/dy))
            
            
            h1term=fun*fun';
            hdxu1=dt*theta*((deriv(1,:)*u1_l+deriv(2,:)*v1_l))*(fun*fun');
            udxh1=dt*theta*(u1I* fun*deriv(1,:)+v1I*fun*deriv(2,:));
            
            
            h0term=h0I*fun;
            a0term=dt*(1-theta)*a0I*fun;
            a1term=dt*theta*a1I*fun;
            hdxu0=-dt*(1-theta)*(deriv(1,:)*u0_l+deriv(2,:)*v0_l)*h0I*fun;
            udxh0=-dt*(1-theta)*(u0I*deriv(1,:)+v0I*deriv(2,:))*h0_l*fun;
            
            km=km+(h1term+hdxu1+udxh1)*weights(Iint)*detJ;
            
            rhs=rhs+(h0term+a0term+a1term+hdxu0+udxh0)*weights(Iint)*detJ;
            
        end % integration points
        
        % assemble global matrix
        
        for i1=1:length(g_l)
            for i2=1:length(g_l)
                kv(g_l(i1),g_l(i2))=kv(g_l(i1),g_l(i2))+km(i1,i2);
            end
            rh(g_l(i1))=rh(g_l(i1))+rhs(i1);
        end
        

    end
   
    disp(' finised assembling matrix ')
     
    [n1L,n2L]=size(L);
    kv=[kv L ; L' zeros(n2L,n2L) ] ; rh=[rh ; Lrhs];
    [n1,n2]=size(kv);
    disp([' number of equations : ',num2str(n1)])
    
    sol=kv\rh;
    
    h1=sol(1:Nnodes) ;
    
    
end   %




