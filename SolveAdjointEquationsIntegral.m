

function [Lambda,Mu,l]=SolveAdjointEquationsIntegral(kv,Lambda,Mu,BoundaryNodes,l,uModel,vModel,uMeas,vMeas,coordinates,connectivity,nip,iteration)
    
    Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
       
    %disp(' Step#2: Solve the adjoint model (right-hand side is a FE integral) ')
    % Step #2: Solve adjoint model
    % -requires change of boundary conditions setting values to zero at boundary
    % -right hand side is different and must be calculated
    
    %disp(' Form right-hand side of ajoint equations ')
    % The right-hand side is:
    dof=2; ndim=2;neqx=Nnodes ;
    neq=dof*Nnodes; rh=zeros(neq,1) ;
    
    [points,weights]=sample('triangle',nip,ndim);
    funInt=cell(nip); derInt=cell(nip);
    for Iint=1:nip
        funInt{Iint}=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        derInt{Iint}=shape_der(Iint,ndim,nod,points);  % dof x nod : dNj/dXi=[dN1/dx dN2/dx dN3/dx; dN1/dy dN2/dy dN3/dy]
    end
    
    for Iele=1:Nele
        con=connectivity(Iele,:);  % nodes of element
        coo=coordinates(con,:) ; % nod x dof =[x1 y1 ; x2 y2 ; x3 y3]
        
        u_l=uModel(connectivity(Iele,:)); v_l=vModel(connectivity(Iele,:)) ;
        uMeas_l=uMeas(connectivity(Iele,:)); vMeas_l=vMeas(connectivity(Iele,:));
        gx_l=connectivity(Iele,:); gy_l=neqx+connectivity(Iele,:);
        b1=zeros(nod,1) ; b2=zeros(nod,1) ;
        
        
        for Iint=1:nip                           % loop over integration points
            fun=funInt{Iint} ; der=derInt{Iint};
            J=der*coo; % (dof x nod) x (nod x dof) = dof x dof
            detJ=det(J);  % det(dof x dof) matrix
            detJw=detJ*weights(Iint);
            % b1=gu * (u-uMeas)
            % b2=gv * (v-vMeas)
            
            
            b1=b1-((u_l-uMeas_l)'*fun)*fun*detJw;
            b2=b2-((v_l-vMeas_l)'*fun)*fun*detJw;
            
        end % integration points
        
        for i1=1:length(gx_l)
            rh(gx_l(i1))=rh(gx_l(i1))+b1(i1);
            rh(gy_l(i1))=rh(gy_l(i1))+b2(i1);
        end
        
    end  % element loop
    
    % modify BC
    %disp(' Define boundary conditions of the adjoint system')
    
    utiedA=[] ; utiedB=[];  vtiedA=[] ; vtiedB=[];
    
    vfixednode=BoundaryNodes;  vfixedvalue=BoundaryNodes*0;
    ufixednode=BoundaryNodes;  ufixedvalue=BoundaryNodes*0;
    [L,Lb]=CreateLuv(Nnodes,ufixednode,ufixedvalue,vfixednode,vfixedvalue,utiedA,utiedB,vtiedA,vtiedB,[],[]);
    
    %  solve: kv [Lambda] = rh    subject to L [Lambda]=Lb
    %            [Mu]                          [Mu]
    % where x0=[Lambda ; Mu] and y=l are starting values if iterative method is used
   
    
    %disp(' Solve adjoint equations '); 
    [sol,l]=solveKApeSymmetric(kv,L,rh,Lb,[Lambda;Mu],l,iteration);
    Lambda=sol(1:Nnodes) ; Mu=sol(Nnodes+1:2*Nnodes);
    
   
    
    