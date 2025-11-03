function TestBoundaryIntegral(coordinates,connectivity,Boundary,h0,u0,v0,u1,v1,a0,a1,dt,CtrlVar)

%%


% The boundary term is:
%
% {h du/dt + u (a-d(hu)/dx-d(hv)/dy)}  n_x N + {h dv/dt + v (a-d(hu)/dx-d(hv)/dy)}  n_y N
%
%

%clc ; clear all
%load TestSAve coordinates connectivity Boundary h etaInt alpha rho rhow g CtrlVar nip
%error('dfas')
% Need to integrate
%  {h du/dt + u (a-d(hu)/dx-d(hv)/dy)}  n_x N_p + {h dv/dt + v (a-d(hu)/dx-d(hv)/dy)}  n_y N_p
%
%  Boundary.ElementsBCu{I} is a list of elements for which Dirichlet is defined along edge I=1:3
%
% loop over each edge
%   loop over elements
%
%    -calculate position s_i and weights of integration points along a 1d line for ndim=1 and nip=nod+1 (at least)
%
%   -determine corresponding points in the 2d (eta,xi) plane, so for example along edge 1
%    well have (s_i,0), along edge 2 it is (0,s_i), along edge 3 it is (1-s ,s )
%
%   -calculate integrand just as done in evaluation the stiffness and the mass matrix
%
%   -sum over integration points with corresponding weights, as usual
%
Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
dof=1; neq=dof*Nnodes; ndim=2; neqx=Nnodes ;


rh=zeros(neq,1);

% find total number of edges in mesh
icount=0;
for iEdge=1:numel(Boundary.Edge)
    icount=icount+numel(union(Boundary.Elements{iEdge},Boundary.Elements{iEdge}));
end

N=nod*nod*icount; Iind=zeros(N,1) ; Jind=Iind ; Xval=Iind;

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

dudt=(u1-u0)/dt ;dvdt=(v1-v0)/dt ;

istak=0;


for iEdge=1:numel(Boundary.Edge)  % loop over edges (just a loop running from 1:3)
    
    [points,weights]=sampleEdge('line',nipEdge,ndim,iEdge);
    % get local coordinates and weights for gamma along the 1d line from -1 to 1
    
    % loop over elements
    for Iele=Boundary.Elements{iEdge}
        
        % gather local quantities from global arrays
        % note the nodal numbering is clockwise!
        con=connectivity(Iele,:);  % nodes of edge of the element
        coo=coordinates(con,:) ; % nod x dof =[x1 y1 ; x2 y2 ; x3 y3]
        
        h0_l=h0(con);  u0_l=u0(con); v0_l=v0(con); u1_l=u1(con); v1_l=v1(con); a0_l=a0(con) ; a1_l=a1(con) ;
        dudt_l=dudt(con); dvdt_l=dvdt(con);
        
        gx_l=con;
        
        c11=zeros(nod,nod) ; b1=zeros(nod,1) ;
        
        for Iint=1:nipEdge                           % loop over integration points
            % the form functions are the ususal 2d form functions
            % but they are evaluated at the integration points along the edge
            fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
            der=shape_der(Iint,ndim,nod,points);  % dof x nod : dNj/dXi=[dN1/dx dN2/dx dN3/dx; dN1/dy dN2/dy dN3/dy]
            J=der*coo; % (dof x nod) x (nod x dof) = dof x dof
            deriv=J\der; % (dof x dof) x (dof x nod) = dof x nod

            
            nxdGamma(1)=der(1,:)*coo(:,2);             %    = J(1,2)                     (1,nod) x (nod,1)=scalar at each integration point
            nxdGamma(2)=-der(2,:)*coo(:,2);            %    =-J(2,2)
            nxdGamma(3)=-(der(1,:)-der(2,:))*coo(:,2); %    =-(J(1,2)-J(2,2))
             
            nydGamma(1)=-der(1,:)*coo(:,1);             % -J(1,1)
            nydGamma(2)=der(2,:)*coo(:,1);              %  J(2,1)
            nydGamma(3)=(der(1,:)-der(2,:))*coo(:,1);   % J(1,1)-J(2,1)
            
            %fprintf('\n \n Iele %i edge %i Iint %i nx %g ny %g \n',Iele,iEdge,Iint,nxdGamma(iEdge),nydGamma(iEdge))
            %connectivity(Iele,:)
            
            u1int=u1_l'*fun; v1int=v1_l'*fun; % scalar
            a1int=a1_l'*fun; a0int=a0_l'*fun; % scalar
            
            du1dx=deriv(1,:)*u1_l;  
            dv1dy=deriv(2,:)*v1_l;  % scalars
            
            
            duintdt=dudt_l'*fun; dvintdt=dvdt_l'*fun;
            
            h0int=h0_l'*fun ; dh0dx=deriv(1,:)*h0_l;  dh0dy=deriv(2,:)*h0_l;  % scalars
            u0int=u0_l'*fun; v0int=v0_l'*fun; % scalar
            du0dx=deriv(1,:)*u0_l;  dv0dy=deriv(2,:)*v0_l;  % scalars
            
                        
            fx=(duintdt*(fun*fun')-u1int*(du1dx*(fun*fun')+u1int*fun*deriv(1,:)+dv1dy*(fun*fun')+v1int*fun*deriv(2,:)))*nxdGamma(iEdge); 
            fy=(dvintdt*(fun*fun')-v1int*(du1dx*(fun*fun')+u1int*fun*deriv(1,:)+dv1dy*(fun*fun')+v1int*fun*deriv(2,:)))*nydGamma(iEdge); 
            c11=c11-(fx+fy)*weights(Iint);
            
            
            fx=(duintdt*h0int-u1int*a1int-u0int*(du0dx*h0int+u0int*dh0dx+dv0dy*h0int+v0int*dh0dy-a0int))*nxdGamma(iEdge)*fun;
            fy=(dvintdt*h0int-v1int*a1int-v0int*(du0dx*h0int+u0int*dh0dx+dv0dy*h0int+v0int*dh0dy-a0int))*nydGamma(iEdge)*fun;
            
            
            b1=b1-(fx+fy)*weights(Iint);
            
            
            
        end % integration points
        
        
        for i1=1:length(gx_l)  ;
            for i2=1:length(gx_l)
                istak=istak+1; Iind(istak)=gx_l(i1); Jind(istak)=gx_l(i2); Xval(istak)=c11(i1,i2);
            end
        end
        
        
        
        for i1=1:length(gx_l)
            rh(gx_l(i1))=rh(gx_l(i1))+b1(i1);
        end
        
    end  % element loop
end


K=sparse(Iind,Jind,Xval,neqx,neqx);

%K=(K+K')/2 ; % I know that the matrix must be symmetric, but numerically this may not be strickly so

%%


end



