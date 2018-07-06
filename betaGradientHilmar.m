function [betaGrad]=betaGradientHilmar(u,v,Lambda,Mu,C,connectivity,coordinates,nip)
    
   % disp(' beta gradient ')
    
    % r_q(u(beta),beta)=0, where
    % r_q= (w_q , L u) = int ( \p x( h eta \p_x u) -\beta^2 u - \rho g h \p_x s) w_q(x) dx   for q=1..N
    % J=I + \lambda_q r_q = I + \lambda_q (w_q , L u)   where u=u_r N_r
    
    % taking derivative with respect to \lambda_i and setting to zero gives (w_i, L u) =0
    
    % taking derivative with respect to u_i of the second term of the augmented const function gives:
    % d/du_i \lambda_q r_q      = d/du_i (\lambda_q (w_q,Lu)= \lambda_q d/du_i (L^T w_q,u_r N_r)
    %                           = \lambda_q d/du_i (L^T w_q,N_i) because I use the Galerking method where
    %                           w_q=N_q and the operation is self adjoint, this is the same system
    
    % taking the derivative with respect to beta_i gives:
    % d/d\beta_i =    \lambda_q   d/d\beta_i (w_q , L u)   where u=u_r N_r
    %            = \lambda_q d/d\beta_i (w_q , \p x( h eta \p_x u) -\beta^2 u - \rho g h \p_x s)
    %            = \lambda_q  (w_q ,  - 2 \beta_r N_r u_r N_r N_i )
    % the inner product is: (w_q ,  - 2 \beta_r N_r u_r N_r N_i ) = -2 int N_q N_i  beta_r N_r u_r N_r , where
    % we sum over r, for each element we have
    % -2 \sum_{p=1}{nip} N_q(x_p) N_I(x_p) \sum_{r=1}^{nod} beta_r(x_p) N_r(x_p) u_r(x_p) N_r(x_p)
    
    Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
    ndim=2; dof=2; neqx=Nnodes ; neqy=Nnodes;
    N=nod*nod*Nele; 
    
    Iindu=zeros(N,1) ; Jindu=Iindu ; Xvalu=Iindu;
    Iindv=zeros(N,1) ; Jindv=Iindv ; Xvalv=Iindv;
    
    [points,weights]=sample('triangle',nip,ndim);
   
    
    
    funInt=cell(nip); derInt=cell(nip);
    for Iint=1:nip
        funInt{Iint}=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        derInt{Iint}=shape_der(Iint,ndim,nod,points);  % dof x nod : dNj/dXi=[dN1/dx dN2/dx dN3/dx; dN1/dy dN2/dy dN3/dy]
    end
    
    istaku=0; istakv=0;
    for Iele=1:Nele
        % gather local quantities from global arrays
        % note the nodal numbering is clockwise!
        con=connectivity(Iele,:);  % nodes of element
        coo=coordinates(con,:) ; % nod x dof =[x1 y1 ; x2 y2 ; x3 y3]
        %h_l=h(con); s_l=s(con); 
        u_l=u(con); v_l=v(con) ; C_l=C(con);
        
        %beta_l= sqrt(C_l.^(-1/m).* (sqrt(u_l.*u_l+v_l.*v_l)).^(1/m-1));
        
        beta_l=sqrt(1./C_l); % only correct for m=1 !!!
        
        gx_l=con; gy_l=neqx+con; c11=zeros(nod,nod) ;c22=c11;

        for Iint=1:nip                           % loop over integration points
            
            fun=funInt{Iint} ; der=derInt{Iint};
            J=der*coo; % (dof x nod) x (nod x dof) = dof x dof
            detJ=det(J);  % det(dof x dof) matrix
            %           deriv=J\der; % (dof x dof) x (dof x nod) = dof x nod
            %           etaI=etaInt(Iele,Iint) ;  % scalar
            %           hI=h_l'*fun ;
            %           CI=C_l'*fun ;
            uI=u_l'*fun ;
            vI=v_l'*fun ;
            betaI=beta_l'*fun;
            
            detJw=detJ*weights(Iint);
            
            c11=c11- 2*betaI*uI*(fun*fun')*detJw;
            c22=c22- 2*betaI*vI*(fun*fun')*detJw;
            
        end % integration points
        
        
        for i1=1:length(gx_l)  ;
            for i2=1:length(gx_l)
                istaku=istaku+1; Iindu(istaku)=gx_l(i1); Jindu(istaku)=gx_l(i2); Xvalu(istaku)=c11(i1,i2);
            end
        end
        
        
        for i1=1:length(gy_l)  ;
            for i2=1:length(gy_l)
                istakv=istakv+1;Iindv(istakv)=gx_l(i1);Jindv(istakv)=gx_l(i2); Xvalv(istakv)=c22(i1,i2);
            end
        end
        
        
    end  % element loop
    
    
    Bu=sparse(Iindu,Jindv,Xvalu,neqx,neqx);
    Bv=sparse(Iindu,Jindv,Xvalv,neqx,neqx);
    
    betaGrad=[Lambda;Mu]'*[Bu;Bv]; betaGrad=betaGrad(:); betaGrad=real(betaGrad);

    betaGrad=EleAverageInterpolate(betaGrad,coordinates,connectivity);
  
    
end




