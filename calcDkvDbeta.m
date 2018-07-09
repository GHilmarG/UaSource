function [DkvDbeta]=calcDkvDbeta(I,s,h,u,v,coordinates,connectivity,Nnodes,Nele,nip,nod,etaInt,gfint,C,m)
    
    ndim=2; dof=2; neq=dof*Nnodes;
    neqx=Nnodes ; neqy=Nnodes;
    
    N=4*nod*nod*Nele; Iind=zeros(N,1) ; Jind=Iind ; Xval=Iind;
    
    [points,weights]=sample('triangle',nip,ndim);
    % get local coordinates and weights

    beta= sqrt(C.^(-1/m).* (sqrt(u.*u+v.*v)).^(1/m-1));
    %x=coordinates(:,1); y=coordinates(:,2); DT = DelaunayTri(x,y); TRI=DT.Triangulation;
    %figure(500) ; trisurf(TRI,x,y,beta) ;  title(' beta ')
    
    dbeta=1e-10;
    beta(I)=beta(I)+sqrt(-1)*dbeta;
    
    funInt=cell(nip); derInt=cell(nip);
    
    for Iint=1:nip
        funInt{Iint}=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        derInt{Iint}=shape_der(Iint,ndim,nod,points);  % dof x nod : dNj/dXi=[dN1/dx dN2/dx dN3/dx; dN1/dy dN2/dy dN3/dy]
    end
    
    istak=0;
    for Iele=1:Nele
        % gather local quantities from global arrays
        % note the nodal numbering is clockwise!
        con=connectivity(Iele,:);  % nodes of element
        coo=coordinates(con,:) ; % nod x dof =[x1 y1 ; x2 y2 ; x3 y3]
        h_l=h(con); s_l=s(con); u_l=u(con); v_l=v(con) ; C_l=C(con);
        beta_l=beta(con);
        
        
        gx_l=con; gy_l=neqx+con;
        
        c11=zeros(nod,nod) ; c12=c11 ; c21=c11 ;c22=c11;
       
        
        for Iint=1:nip                           % loop over integration points
            
            
            
            fun=funInt{Iint} ; der=derInt{Iint};
            
            J=der*coo; % (dof x nod) x (nod x dof) = dof x dof
            detJ=det(J);  % det(dof x dof) matrix
            deriv=J\der; % (dof x dof) x (dof x nod) = dof x nod
            etaI=etaInt(Iele,Iint) ;  % scalar
            hI=h_l'*fun ;
            CI=C_l'*fun ;
            uI=u_l'*fun ;
            vI=v_l'*fun ;
            betaI=beta_l'*fun;
            
            
            
            detJw=detJ*weights(Iint);
            
            d1d1=deriv(1,:)'*deriv(1,:);
            d2d2=deriv(2,:)'*deriv(2,:);
            c11=c11-(etaI*hI*(4*d1d1+d2d2)+ betaI^2*(fun*fun'))*detJw;
            c12=c12-(etaI*hI*(2*deriv(1,:)'*deriv(2,:)+deriv(2,:)'*deriv(1,:)))*detJw;
            %c21=c21-(etaI*hI*(2*deriv(2,:)'*deriv(1,:)+deriv(1,:)'*deriv(2,:)))*detJw;
            c21=c12;
            c22=c22-(etaI*hI*(4*d2d2+d1d1)+betaI^2*(fun*fun'))*detJw;
            
     
        end % integration points
        
        for i1=1:length(gx_l)  ;
            for i2=1:length(gx_l)
                istak=istak+1; Iind(istak)=gx_l(i1); Jind(istak)=gx_l(i2); Xval(istak)=c11(i1,i2);
            end
        end
        
        
        for i1=1:length(gy_l)  ;
            for i2=1:length(gy_l)
                istak=istak+1;
                Iind(istak)=gy_l(i1);
                Jind(istak)=gy_l(i2);
                Xval(istak)=c22(i1,i2);
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
        
    
    end  % element loop
    
    DkvDbeta=sparse(Iind,Jind,Xval,neqx+neqy,neqx+neqy);
    %DkvDbeta=(DkvDbeta+DkvDbeta.')/2 ;
    
    DkvDbeta=imag(DkvDbeta)/dbeta;
    
end



