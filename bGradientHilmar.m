function [bGrad]=bGradientHilmar(u,v,s,h,etaInt,gfint,rho,rhow,g,alpha,m,C,Lambda,Mu,connectivity,coordinates,nip)
    
     % calculates \lambda_q \p r_q/\p I =Lambda b_i \frac{\p}{\b_i}  <L u-f| N_q> 
    % where the FE formulation used does NOT have  Weertman BC and the natural boundary condition
    % Note: this has not been vectorized
    
    
    % disp(' b gradient ')
    TestGradient=0 ; % tests gradient by calculating it numerically as well
    ReturnTestGradient=0;
    
    
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
        u_l=u(con); v_l=v(con) ; s_l=s(con);
        
        
        gx_l=con; gy_l=neqx+con; Fx=zeros(nod,nod) ;Fy=Fx;
        
        for Iint=1:nip                           % loop over integration points
            
            fun=funInt{Iint} ; der=derInt{Iint};
            J=der*coo; % (dof x nod) x (nod x dof) = dof x dof
            detJ=det(J);  % det(dof x dof) matrix
            deriv=J\der; % (dof x dof) x (dof x nod) = dof x nod
            etaI=etaInt(Iele,Iint) ;  % scalar
            %           hI=h_l'*fun ;
            %           CI=C_l'*fun ;
            %uI=u_l'*fun ;
            %vI=v_l'*fun ;
            
            dudxI=deriv(1,:)*u_l;
            dudyI=deriv(2,:)*u_l;
            dvdxI=deriv(1,:)*v_l;
            dvdyI=deriv(2,:)*v_l;
            
            
            detJw=detJ*weights(Iint);
            
                        
            Fx=Fx+((4*dudxI+2*dvdyI)*etaI*deriv(1,:)'*fun'+(dvdxI+dudyI)*etaI*deriv(2,:)'*fun'+...
                rho*g*(deriv(1,:)*s_l*cos(alpha)-sin(alpha))*(fun*fun'))*detJw;
            
            Fy=Fy+((4*dvdyI+2*dudxI)*etaI*deriv(2,:)'*fun'+(dudyI+dvdxI)*etaI*deriv(1,:)'*fun'+...
                rho*g*deriv(2,:)*s_l*cos(alpha)*(fun*fun'))*detJw;
            
            
            
            
            
            
        end % integration points
        
        
        for i1=1:length(gx_l)  ;
            for i2=1:length(gx_l)
                istaku=istaku+1; Iindu(istaku)=gx_l(i1); Jindu(istaku)=gx_l(i2); Xvalu(istaku)=Fx(i1,i2);
            end
        end
        
        
        for i1=1:length(gx_l)  ; % I must use gx rather than gy because I'm putting this into two seperate vectors
            for i2=1:length(gx_l)
                istakv=istakv+1;Iindv(istakv)=gx_l(i1);Jindv(istakv)=gx_l(i2); Xvalv(istakv)=Fy(i1,i2);
            end
        end
        
        
        
    end  % element loop
    
    
    Bu=sparse(Iindu,Jindv,Xvalu,neqx,neqx);
    Bv=sparse(Iindu,Jindv,Xvalv,neqx,neqx);
    
    bGrad=[Lambda;Mu]'*[Bu';Bv']; bGrad=bGrad(:); bGrad=real(bGrad);
    
    bGrad=EleAverageInterpolate(bGrad,coordinates,connectivity);
    
    if TestGradient==1
        disp(' testing d kv / db  ' )  % Takes a lot of time, only good for verification purposes
        bGradNumeric=zeros(Nnodes,1);
        bGradNumeric1=zeros(Nnodes,1);
        bGradNumeric2=zeros(Nnodes,1);
        for I=1:Nnodes
            
            if mod(10*I,Nnodes)<10 ; disp([num2str(100*I/Nnodes),' % ']) ; end
            [DkvDb,DrhDb]=calcDkvDbComplex(I,s,h,u,v,coordinates,connectivity,Nnodes,Nele,nip,nod,etaInt,gfint,C,m,alpha,rho,rhow,g);
            
            bGradNumeric(I)=-[Lambda;Mu]'*DkvDb*[u;v]+[Lambda;Mu]'*DrhDb;
            bGradNumeric1(I)=-[Lambda;Mu]'*DkvDb*[u;v];
            bGradNumeric2(I)=[Lambda;Mu]'*DrhDb;
        end
        
        bGradNumeric=EleAverageInterpolate(bGradNumeric,coordinates,connectivity);
        bGradNumeric1=EleAverageInterpolate(bGradNumeric1,coordinates,connectivity);
        bGradNumeric2=EleAverageInterpolate(bGradNumeric2,coordinates,connectivity);
       
      
        
        
        x=coordinates(:,1); y=coordinates(:,2); DT = DelaunayTri(x,y); TRI=DT.Triangulation;
        figure(600) ; trisurf(TRI,x,y,bGradNumeric) ;  title(' bGradNumeric ')
        figure(601) ; trisurf(TRI,x,y,bGrad) ;  title(' bGrad ')
        figure(602) ; trisurf(TRI,x,y,bGradNumeric-bGrad) ;  title(' bGradNumeric-bGrad ')
        figure(603) ; trisurf(TRI,x,y,bGradNumeric1) ;  title(' bGradNumeric1 ')
        figure(604) ; trisurf(TRI,x,y,bGradNumeric2) ;  title(' bGradNumeric2 ')
        
        figure(605) ; trisurf(TRI,x,y,bGradNumeric1-bGrad) ;  title(' bGradNumeric1-bGrad ')
        figure(606) ; trisurf(TRI,x,y,bGradNumeric2-bGrad) ;  title(' bGradNumeric2-bGrad ')
        
        save Testdb TRI x y bGrad bGradNumeric
        
        if ReturnTestGradient==1
            bGrad=bGradNumeric;
        end
        
    end
    
    
    
end




