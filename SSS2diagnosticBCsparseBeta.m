function [u,v,exx,eyy,exy,etaInt,beta2,xint,yint,kv,lambdau]=...
        SSS2diagnosticBCsparseBeta(s,h,u,v,coordinates,connectivity,Nnodes,Nele,nip,nod,etaInt,AGlen,C,...
        L,Lb,lambdau,n,m,alpha,rho,rhow,g,Itime)
    
   
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
    disp([' taub = ',num2str(taub),' u = ',num2str(mean(C)*taub^m)])
    
    [points,weights]=sample('triangle',nip,ndim);
    % get local coordinates and weights
    
    diff=1e10; 
    
    %[--  constanst related to non-linear loop. Values used might affect solution 
    beta2Imax=1e5; tol=1e-5 ;  etamax=1e20 ; % 
    %]-
    
    % beta2= C^(-1/m)* (u.*u+v.*v)^(1/m-1) ;
    
    if n==1 && m==1 ; 
        iteration_max=1 ; iteration_min=1 ;
    else
        disp(' non-linear ')
        iteration_max=50; iteration_min=2 ;
    end
    
    iteration=0;
    
    funInt=cell(1,3); derInt=cell(1,3);
    for Iint=1:nip
        funInt{Iint}=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        derInt{Iint}=shape_der(Iint,ndim,nod,points);  % dof x nod : dNj/dXi=[dN1/dx dN2/dx dN3/dx; dN1/dy dN2/dy dN3/dy]
    end
        
    while (diff > tol  && iteration < iteration_max )  || iteration < iteration_min
        iteration=iteration+1;
        disp([' SSS iteration # : ',num2str(iteration)])
        rh=zeros(neq,1) ; 
        
        istak=0;
        % element loop
        tStartAssembly=tic;
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
                beta2(Iele,Iint)=beta2I;
                
                % if beta2I > beta2Imax ; beta2I=beta2Imax ; end
                
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
                
                detJw=detJ*weights(Iint);
                d1d1=deriv(1,:)'*deriv(1,:);
                d2d2=deriv(2,:)'*deriv(2,:);
                c11=c11-(etaI*hI*(4*d1d1+d2d2)+ beta2I*(fun*fun'))*detJw;
                c12=c12-(etaI*hI*(2*deriv(1,:)'*deriv(2,:)+deriv(2,:)'*deriv(1,:)))*detJw;
                %c21=c21-(etaI*hI*(2*deriv(2,:)'*deriv(1,:)+deriv(1,:)'*deriv(2,:)))*detJw;
                c21=c12;
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
            
            
            for i1=1:length(gx_l)
                rh(gx_l(i1))=rh(gx_l(i1))+b1(i1);
                rh(gy_l(i1))=rh(gy_l(i1))+b2(i1);
            end
           
        end  % element loop

       
        
        kv=sparse(Iind,Jind,Xval,neqx+neqy,neqx+neqy);
    
        kv=(kv+kv')/2 ; % I know that the matrix must be symmetric, but numerically this may not be strickly so
        
         tElapsedAssembly=toc(tStartAssembly);
        disp([' time for matrix assembly ',num2str(tElapsedAssembly)])
        
        
        tStartSolve=tic; 
        [sol,lambdau]=solveKApeSymmetric(kv,L,rh,Lb,[u;v],lambdau,iteration+Itime-1);
        tElapsedSolve=toc(tStartSolve);
        disp([' time for matrix solution ',num2str(tElapsedSolve)])
        
        ulast=u ; vlast=v ;
        u=sol(1:Nnodes) ; v=sol(Nnodes+1:2*Nnodes);
        
       
        
        
        diff=(max(abs(u-ulast))+max(abs(v-vlast)))/(max(abs(u))+max(abs(v)));
        if iteration>1 ;disp(['diff : ',num2str(diff)]) ; end
     
        
     
          
        if calcStrainRates==1
             tStartStrainRates=tic; 
            %disp(' start calculating strain rates and the effective strain rate')
            %[---   get strain rates and calculate effective viscosity
            exx=zeros(Nele,nip); eyy=exx ;exy=exx;
            xint=zeros(Nele,nip); yint=xint;
            
            for Iele=1:Nele                           % loop over elements
                
                con=connectivity(Iele,:);  % nodes of element
                coo=coordinates(con,:) ; % nod x dof =[x1 y1 ; x2 y2 ; x3 y3]
                u_l=u(con); v_l=v(con) ; 

                for Iint=1:nip                           % loop over integration points
                    
                    %fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
                    %der=shape_der(Iint,ndim,nod,points);  % dof x nod : dNj/dXi=[dN1/dx dN2/dx dN3/dx; dN1/dy dN2/dy dN3/dy]
                    fun=funInt{Iint} ; der=derInt{Iint};
                    
                    J=der*coo; % (dof x nod) x (nod x dof) = dof x dof
                    %detJ=det(J); 
                    %iJ=inv(J); % dof x dof matrix
                    deriv=J\der; % (dof x dof) x (dof x nod) = dof x nod
                    %hI=h_l'*fun ;  CI=C_l'*fun ; uI=u_l'*fun ; vI=v_l'*fun ;J=der*coord;
                   
                    
                                       
                    exx(Iele,Iint)=deriv(1,:)*u_l;       
                    eyy(Iele,Iint)=deriv(2,:)*v_l;       
                    exy(Iele,Iint)=0.5*(deriv(2,:)*u_l+deriv(1,:)*v_l);
                    xint(Iele,Iint)=coo(:,1)'*fun ; yint(Iele,Iint)=coo(:,2)'*fun ;
                    
                    
                    
                    e=sqrt(exx(Iele,Iint)^2+eyy(Iele,Iint)^2+exx(Iele,Iint)*eyy(Iele,Iint)+exy(Iele,Iint)^2); % consider taking out of loop
                    etaInt(Iele,Iint)=0.5*AGlen^(-1/n)*e^((1-n)/n);
                    
                    
                end
            end
            etaInt(etaInt>etamax)=etamax;
            %disp(' finished calculating strain rates and the effective strain rate')
            tElapsedStrainRates=toc(tStartStrainRates);
            disp([' time for strain rate calculation  ',num2str(tElapsedStrainRates)])
            
        end
       
    end   
    
end



