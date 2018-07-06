function [u,v]=SSS2diagnostic(s,h,u,v,coordinates,connectivity,Nnodes,Nele,nip,nod,etaInt,AGlen,C,uzero,vzero,L,Lb,n,m,alpha,rho,g)
              
    'fdasdfa'
    [n1L,n2L]=size(L);
    
    
    ndim=2; dof=2; neq=dof*Nnodes;
    neqx=Nnodes ; neqy=Nnodes;
    
    
    % solves the non-linear 2d shallow ice stream equation
    % and shallow shelf equation
    % using the method of finite elements
    %
    
    rhog=rho*g; ca=cos(alpha); sa=sin(alpha);
    taub=rhog*sin(alpha)*mean(h);
    disp([' taub = ',num2str(taub),' u = ',num2str(mean(C)*taub^m)])
    
    [points,weights]=sample('triangle',nip,ndim);
    % get local coordinates and weights
    
    diff=1e10; 
    
    %[--  constanst related to non-linear loop. Values used might affect solution 
    beta2Imax=1e5; tol=1e-3 ;  % 
    %]-
    
    % beta2= C^(-1/m)* (u.*u+v.*v)^(1/m-1) ;
    
    if n==1 && m==1 ; 
        iteration_max=1 ; iteration_min=2 ;
    else
        disp(' non-linear ')
        iteration_max=50; iteration_min=2 ;
    end
    
    iteration=0;
    
    while (diff > tol  && iteration < iteration_max )  || iteration < iteration_min
        iteration=iteration+1;
        disp([' SSS iteration # : ',num2str(iteration)])
        rh=zeros(neq,1) ; kv=zeros(neq,neq)  ;
        
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
                
                fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
                der=shape_der(Iint,ndim,nod,points);  % dof x nod : dNj/dXi=[dN1/dx dN2/dx dN3/dx; dN1/dy dN2/dy dN3/dy]
                
                J=der*coo; % (dof x nod) x (nod x dof) = dof x dof
                detJ=det(J); iJ=inv(J); % dof x dof matrix
                deriv=iJ*der; % (dof x dof) x (dof x nod) = dof x nod
                etaI=etaInt(Iele,Iint) ;  % scalar
                hI=h_l'*fun ;
                CI=C_l'*fun ;
                uI=u_l'*fun ;
                vI=v_l'*fun ;
                
                beta2I= CI^(-1/m)* (sqrt(uI.*uI+vI.*vI))^(1/m-1) ; % scalar
                
                
                if m~=1 ; beta2I=min([beta2I ; beta2Imax]) ; end
               
                
                % h_l , s_l , beta_l : nod x 1
                % h_l' * fun
                % hI=h_i N_i = h_l' * fun
                
                % X component:
                % left-hand side
                % -(4 h etaI u_x + 2 h etaI v_y ) w_x - h etaI (v_x + u_y ) w_y
                % -( 4 hI etaI u_p dNp/dx dNq/dx +  2 h etaI v_p dNp/dy dNq/dx + hI etaI (v_p dNp/dx+  u_p dNp/dy ) dNq/dy
                %= -etaI hI (4 dNp/dx dNq/dx + dNp/dy dNq/dy) u_p  - etaI hI (2 dNp/dy dNq/dx + dNp/dx dNq/dy) v_p
                %=   -etaI hI (4 deriv(1,:)'*deriv(1,:) + deriv(2,:)'*deriv(2,:) ) u_p
                %    -etaI hI (2 deriv(1,:)'*deriv(2,:) + deriv(2,:)'*deriv(1,:) ) v_p
                %
                % tbx =beta2 u
                % where beta2= C^(-1/(2*m))* (u_l.*u_l+v_l.*v_l)^(1/m-1) is nodx1
                % tbx w = beta2'*fun * N_q N_p u_p = (beff'*fun) * (fun*fun') u_p
                %
                
                
                % right hand side
                % rho g h (p_x s \cos \alpha - \sin \alpha) w
                % rho g hI (s_q dNq/dx cos(alpha) - sin(alpha)) N_p
                % rho g hI (deriv(1,:)*s_l cos(alpha) - sin(alpha))* fun
                
                % Y component:
                % left-hand side
                % -(4 h etaI v_y + 2 h etaI u_x ) w_y - h etaI (u_y + v_x ) w_x
                % -( 4 hI etaI v_p dNp/dy dNq/dy +  2 h etaI u_p dNp/dx dNq/dy + hI etaI (u_p dNp/dy+  v_p dNp/dx ) dNq/dx
                %= -etaI hI (4 dNp/dy dNq/dy + dNp/dx dNq/dx) v_p  - etaI hI (2 dNp/dx dNq/dy + dNp/dy dNq/dx) u_p
                %=   -etaI hI (4 deriv(2,:)'*deriv(2,:) + deriv(1,:)'*deriv(1,:) ) v_p
                %    -etaI hI (2 deriv(2,:)'*deriv(1,:) + deriv(1,:)'*deriv(2,:) ) u_p
                %
                % tby =beta2 v
                % where beta2= C^(-1/(2*m))* (u_l.*u_l+v_l.*v_l)^(1/m-1) is nodx1
                % tby w = beta2'*fun * N_q N_p v_p = (beta2'*fun) * (fun*fun') v_p
                
                %
                % right hand side
                % rho g hI p_y s \cos \alpha  w
                % rho g hI s_q dNq/dy cos(alpha) N_p
                % rho g hI deriv(2,:)*s_l*cos(alpha)* fun
                
                
                % X-LHS summary:
                %=   -etaI hI (4 deriv(1,:)'*deriv(1,:) + deriv(2,:)'*deriv(2,:) ) u_p
                %    -etaI hI (2 deriv(1,:)'*deriv(2,:) + deriv(2,:)'*deriv(1,:) ) v_p
                %   - beta2I * (fun*fun') u_p
                % X-RHS summary:
                % rho g hI (deriv(1,:)*s_l cos(alpha) - sin(alpha))* fun
                
                % Y-LHS summary:
                %=   -etaI hI (4 deriv(2,:)'*deriv(2,:) + deriv(1,:)'*deriv(1,:) ) v_p
                %    -etaI hI (2 deriv(2,:)'*deriv(1,:) + deriv(1,:)'*deriv(2,:) ) u_p
                %   - beta2I * (fun*fun') v_p
                % Y-RHS summary:
                % rho g hI deriv(2,:)*s_l cos(alpha) * fun
                
                % [c11 c12] [u] = [b1]
                % [c21 c22] [v]   [b2]
                % with
                % c11=  -etaI hI(4 deriv(1,:)'*deriv(1,:) + deriv(2,:)'*deriv(2,:)) - beta2I * (fun*fun')
                % c12=  -etaI hI (2 deriv(1,:)'*deriv(2,:) + deriv(2,:)'*deriv(1,:) )
                % c21=  -etaI hI (2 deriv(2,:)'*deriv(1,:) + deriv(1,:)'*deriv(2,:) )
                % c22=  -etaI hI (4 deriv(2,:)'*deriv(2,:) + deriv(1,:)'*deriv(1,:) ) - beta2I * (fun*fun')
                %
                % b1= rho g hI (deriv(1,:)*s_l cos(alpha) - sin(alpha))* fun
                % b2= rho g hI deriv(2,:)*s_l cos(alpha) * fun
                
                c11=c11-(etaI*hI*(4*deriv(1,:)'*deriv(1,:)+deriv(2,:)'*deriv(2,:))+ beta2I*(fun*fun'))*detJ*weights(Iint);
                c12=c12-(etaI*hI*(2*deriv(1,:)'*deriv(2,:)+deriv(2,:)'*deriv(1,:)))*detJ*weights(Iint);
                c21=c21-(etaI*hI*(2*deriv(2,:)'*deriv(1,:)+deriv(1,:)'*deriv(2,:)))*detJ*weights(Iint);
                c22=c22-(etaI*hI*(4*deriv(2,:)'*deriv(2,:)+deriv(1,:)'*deriv(1,:))+beta2I*(fun*fun'))*detJ*weights(Iint);
                       
                b1=b1+rhog*hI*((deriv(1,:)*s_l)*ca-sa)*fun*detJ*weights(Iint);
                b2=b2+rhog*hI*(deriv(2,:)*s_l)*ca*fun*detJ*weights(Iint);
        
        
            end % integration points
            
            %c11
            % assemble global matrix
            for i1=1:length(gx_l)
                for i2=1:length(gx_l)
                    kv(gx_l(i1),gx_l(i2))=kv(gx_l(i1),gx_l(i2))+c11(i1,i2);
                end
            end
            for i1=1:length(gy_l)
                for i2=1:length(gy_l)
                    kv(gy_l(i1),gy_l(i2))=kv(gy_l(i1),gy_l(i2))+c22(i1,i2);
                end
            end
            for i1=1:length(gx_l)
                for i2=1:length(gy_l)
                    kv(gx_l(i1),gy_l(i2))=kv(gx_l(i1),gy_l(i2))+c12(i1,i2);
                end
            end
            for i1=1:length(gy_l)
                for i2=1:length(gx_l)
                    kv(gy_l(i1),gx_l(i2))=kv(gy_l(i1),gx_l(i2))+c21(i1,i2);
                end
            end
            for i1=1:length(gx_l)
                rh(gx_l(i1))=rh(gx_l(i1))+b1(i1);
                rh(gy_l(i1))=rh(gy_l(i1))+b2(i1);
            end
        end  % element loop
         
        % BC
        % use penalty to set nodes uzero and vzero to zero
        pen=1e20;
        if ~isempty(uzero)
            for II=1:length(uzero)
                kv(uzero(II),uzero(II))=pen;
            end
        end
        
        if ~isempty(vzero)
            for II=1:length(vzero)
                kv(vzero(II)+neqx,vzero(II)+neqx)=pen;
            end
        end
        
        
        
        kv=[kv L ; L' zeros(n2L,n2L) ] ; rh=[rh ; Lb];
        [n1,n2]=size(kv);
        disp([' number of equations : ',num2str(n1)])
        sol=sparse(kv\rh);
        ulast=u ; vlast=v ;
        u=sol(1:Nnodes) ; v=sol(Nnodes+1:2*Nnodes);
        
        diff=(max(abs(u-ulast))+max(abs(v-vlast)))/(max(abs(u))+max(abs(v)));
        disp(['diff : ',num2str(diff)])
        
    end   %
    
end
    
    
    
