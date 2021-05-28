function [Deriv,detJ]=derivVector(coordinates,connectivity,nip,points,Iint)

narginchk(5,5)

% calculates the derivatives of form functions with respect to x and y and the Jacobian
% at a given integration point (This function must be called within a loop over integration points)

% Deriv : Nele x dof x nod
%  detJ : Nele

%     To calculate the x y derivatives of f where f is defined on nodes, at integration points do:
%
%   fnod=reshape(f(connectivity,1), exx=zeros(Nele,nip);
%   coox=reshape(coordinates(connectivity,1),Nele,nod);
%   cooy=reshape(coordinates(connectivity,2),Nele,nod);
%   dfx=zeros(Nele,nod); dfy=zeros(Nele,nod);
%   for Iint=1:nip                           % loop over integration points
%         fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
%         xint(:,Iint)=coox*fun;
%         yint(:,Iint)=cooy*fun;
%         fint(:,Iint)=fnod*fun;
%         [Deriv]=derivVector(coordinates,connectivity,nip,Iint); % Nele x dof x nod
%         for I=1:nod
%             dfx(:,Iint)=dfx(:,Iint)+Deriv(:,1,I).*fnod(:,I);
%             dfy(:,Iint)=dfy(:,Iint)+Deriv(:,2,I).*fnod(:,I);
%         end
%   end
%
% A simple way of calculating the derivatives the nodal variable f is using:
% [dfdx,dfdy,xint,yint]=calcFEderivatives(f,coordinates,connectivity,nip,CtrlVar)
%



[Nele,nod]=size(connectivity); ndim=2; dof=2;
%[points,weights]=sample('triangle',nip,ndim);


if Nele==0
    Deriv=[];
    detJ=[];
    return
end

%hnod=reshape(h(connectivity,1),Nele,nod);
coox=reshape(coordinates(connectivity,1),Nele,nod);
cooy=reshape(coordinates(connectivity,2),Nele,nod);


%hint=zeros(Nele,nip) ;
%J11=zeros(Nele) ; J12=zeros(Nele) ; J21=zeros(Nele) ; J22=zeros(Nele);
%iJ11=zeros(Nele) ; iJ12=zeros(Nele) ; iJ21=zeros(Nele) ; iJ22=zeros(Nele);



% fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
% N=N(r,s)
der=shape_der(Iint,ndim,nod,points);  % dof x nod : dNj/dXi=[dN1/dr dN2/dr dN3/dr; dN1/ds dN2/ds dN3/ds]
%hint=hnod*fun;

% calculate J, detJ, and deriv in vectorized way
%  x=x_p N_p   ; h=h_q N_q
% der=dN_p/dxi_q
% J11=dNp/dxi_1 x_p
% (x,y) and N(r,s)
% x=x_p N_p(r,s) , y=y_p N_p(r,s) , h=h_p N_p(r,s)
%
%
% dh/dx=dh/dr dr/dx+ dh/ds ds/dx =  h_p (dN_p/dr dr/dx + dN_p/ds ds/dx)
% dh/dy=dh/dr dr/dy+ dh/ds ds/dy =  h_p (dN_p/dr dr/dy + dN_p/ds ds/dy)
%
% or
%
% [dh/dx]= h_p [dr/dx   ds/dx ] [dN_p/dr]
% [dh/dy]      [dr/dy   ds/dy ] [dN_p/ds]
%        = h_p J^{-1}  [dN_p/dr ; dN_p/s]
%        =                    J^{-1} [ dN1/dr  dN2/dr ... dNnod/dr ] [h_1]
%                                    [ dN1/ds  dN2/ds ... dNnod/ds ] [h_2]
%                                                                    [   ]
%                                                                    [hnod]
%        = J^{-1} der h
%        =deriv h
%
%  der(1:dof,1:nod)=[ dN1/dr  dN2/dr ... dNnod/dr ]
%                   [ dN1/ds  dN2/ds ... dNnod/ds ]
%
% The Jakobian is: J(1:dof,1:dof)= [dx/dr  dy/dr ]
%                                  [dx/ds  dy/ds ]
% and can be calculated as
% J=der*coo; % (dof x nod) x (nod x dof) = dof x dof
% this expression is only evaluated at the integration points,
% der is the derivative of form functions evaluated at integration points
%
% At each integration point
% J11(1:Nele)=dx/dr=x_p(1:Nele) dN_p/dr = coox*der(1,:)'
% where coox(1:Nele,1:nod)
% J12(1:Nele)=dy/dr=cooy*der(1,:)'
% deriv=inv(J)*der=iJ(1:dof,1:dof) der(1:dof,1:nod) = dof x nod
%
% deriv=[ dN1/dx dN2/dx ... dNnod/dx ] = [dr/dx ds/dx] [dN1/dr   dN2/dr ... dNnod/dr ]
%       [ dN1/dy dN2/dy ... dNnod/dy ]   [dr/dy ds/dy] [dN1/ds   dN2/ds ... dNnod/ds ]
%
% iJac=[dr/dx ds/dx]
%      [dr/dy ds/dy]
%
% iJac=[J22 -J12]  / detJ
%      [-J21 J22]
%
% dN1/dx=iJ11*der(1,1)+iJ12*der(2,1)
% dN2/dx=iJ11*der(1,2)+iJ12*der(2,2)
% dNnod/dx=iJ11*der(1,nod)+iJ12*der(2,nod)
% dNnod/dy=iJ21*der(2,nod)+iJ22*der(2,nod)
% deriv=iJ der; % (dof x dof) x (dof x nod) = dof x nod
%

% this is the Jacobian
J11=coox*der(1,:)' ; % dN1/dxi1 : (Nele x nod) x nod = Nele
J12=cooy*der(1,:)' ;
J21=coox*der(2,:)' ;
J22=cooy*der(2,:)' ;
detJ=J11.*J22-J12.*J21; % This is the determinant of the Jacobian

Deriv=zeros(Nele,dof,nod);

for J=1:nod
    
    Deriv(:,1,J)=(J22*der(1,J)-J12*der(2,J))./detJ;
    Deriv(:,2,J)=(-J21*der(1,J)+J11*der(2,J))./detJ;
    
end

% clearvars J11 J12 J21 J22
% 
% if all(detJ<0) 
%     warning('derivVector:AllElementsInsideOut','Negative determinant in all elements \n  ') ;
% elseif any(detJ<0) 
%     warning('derivVector:SomeElementsInsideOut','Negative determinant in some elements \n  ') ;
% end

end




