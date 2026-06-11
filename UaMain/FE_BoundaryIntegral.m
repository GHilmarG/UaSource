function [K,rhs]=FE_BoundaryIntegral(CtrlVar,coordinates,connectivity,fxx,fyy,fxy,fyx,bx,by)


% do boundary integral <fxx ,N_q> u_p + <fxy, N_q>  v_q = <bx,N_q>

% just testing



Nodes=max(connectivity(:)); [Nele,nod]=size(connectivity);


neq=2*Nodes; neqx=Nodes; ndim=2;

rhs=zeros(neq,1);

nEdges=3; % total number of edges in mesh
N=nod*nod*nEdges;

I=4*zeros(N,1) ; J=4*zeros(N,1);

%Ixy=zeros(N,1) ; Jxy=zeros(N,1);
fval=4*zeros(N,1);
counter=0;

switch CtrlVar.TriNodes
    case 3 % 1 exact for linear variatoin
        nipEdge=1;
    case 6   % 3 exact for second degree polynomials
        nipEdge=4;
    case 10 % mini
        nipEdge=4;
    otherwise
        error(' case not recognised, TriNodes value incorrect')
end


% loop over elements
for Iele=1:Nele
    
    % gather local quantities from global arrays
    % note the nodal numbering is clockwise!
    con=connectivity(Iele,:);  % nodes of edge of the element
    coo=coordinates(con,:) ; % nod x dof =[x1 y1 ; x2 y2 ; x3 y3]
    
    fxx_l=fxx(con) ; fyy_l=fyy(con) ;
    fxy_l=fxy(con) ; fyx_l=fyx(con) ;
    bx_l=bx(con) ; by_l=by(con) ;
    
    gx_l=con; gy_l=con+neqx;
    
    Kxx=zeros(nod,nod) ;  Kyy=zeros(nod,nod) ;
    Kxy=zeros(nod,nod) ;  Kyx=zeros(nod,nod) ;
    Bx=zeros(nod,1) ;    By=zeros(nod,1) ;
    
    
    for iEdge=1:3  % loop over edges (just a loop running from 1:3)
        
        [points,weights]=sampleEdge('line',nipEdge,ndim,iEdge);
        % get local coordinates and weights for gamma along the 1d line from 0 to 1
        
        for Iint=1:nipEdge                           % loop over integration points
            % the form functions are the ususal 2d form functions
            % but they are evaluated at the integration points along the edge
            fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
            der=shape_der(Iint,ndim,nod,points);  % dof x nod : dNj/dXi=[dN1/dx dN2/dx dN3/dx; dN1/dy dN2/dy dN3/dy]
            Jac=der*coo; % (dof x nod) x (nod x dof) = dof x dof
            deriv=Jac\der; % (dof x dof) x (dof x nod) = dof x nod
            
            
            nxdGamma(1)=der(1,:)*coo(:,2); % (1,nod) x (nod,1)=scalar at each integration point
            nxdGamma(2)=-der(2,:)*coo(:,2);
            nxdGamma(3)=-(der(1,:)-der(2,:))*coo(:,2);
            
            nydGamma(1)=-der(1,:)*coo(:,1);
            nydGamma(2)=der(2,:)*coo(:,1);
            nydGamma(3)=(der(1,:)-der(2,:))*coo(:,1);
            
            % [nxdGamma ;nydGamma]
            
            
            
            % sum over nodes
            fxxint=fxx_l'*fun ;   fxyint=fxy_l'*fun ;
            fyyint=fyy_l'*fun ;   fyxint=fyx_l'*fun ;
            
            bxint=bx_l'*fun; byint=by_l'*fun;
            
            Fxx=fxxint*(fun*fun')*nxdGamma(iEdge); Fxy=fxyint*(fun*fun')*nxdGamma(iEdge);
            Fyy=fyyint*(fun*fun')*nydGamma(iEdge); Fyx=fyxint*(fun*fun')*nydGamma(iEdge);
            
            Kxx=Kxx+Fxx*weights(Iint);  Kxy=Kxy+Fxy*weights(Iint);
            Kyy=Kyy+Fyy*weights(Iint);  Kyx=Kyx+Fyx*weights(Iint);
            
            
            
            rhx=bxint*nxdGamma(iEdge)*fun;
            rhy=byint*nydGamma(iEdge)*fun;
            
            Bx=Bx+rhx*weights(Iint);
            By=By+rhy*weights(Iint);
            %
            
        end % integration points
        
    end
    for i1=1:length(gx_l)  ;
        for i2=1:length(gx_l)
            counter=counter+1; I(counter)=gx_l(i1); J(counter)=gx_l(i2); fval(counter)=Kxx(i1,i2);
            counter=counter+1; I(counter)=gx_l(i1); J(counter)=gy_l(i2); fval(counter)=Kxy(i1,i2);
            counter=counter+1; I(counter)=gy_l(i1); J(counter)=gx_l(i2); fval(counter)=Kyx(i1,i2);
            counter=counter+1; I(counter)=gy_l(i1); J(counter)=gy_l(i2); fval(counter)=Kyy(i1,i2);
        end
    end
    
    for i1=1:length(gx_l)
        rhs(gx_l(i1))=rhs(gx_l(i1))+Bx(i1);
        rhs(gy_l(i1))=rhs(gy_l(i1))+By(i1);
    end
    
    
end  % element loop

K=sparse(I,J,fval,neq,neq);
end


