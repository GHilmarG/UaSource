function [ud,vd]=uvSSHEETplus(CtrlVar,MUA,BCs,AGlen,n,h,txzb,tyzb)

%  calculates deformational velocity based on the SSHEET (SIA) approximation
%  
%  u and v are nodal velocities
%
%  u= -(2A/(n+1)) (rho g)^n | grad_{xy} s|^(n-1) h^(n+1) \p_x s
%
% u= (2A/(n+1))   (txzb^2+tyxz^2)^((n-1)/2)  \txzb  h 
%
%  N_p N_q u_q = -N_p (2A/(n+1)) (rho g)^n | grad_{xy} s|^(n-1) h^(n+1) \p_x s
%
% This is a linear system and I just need to solve
%
% [M Lu'] [u]  = [rhs  ]
% [Lu 0 ] [lu] = [Lrhsu]
%
% All matrices are Nnodes x Nnodes, apart from:
% Luv is #uv constraints x 2 Nnodes
%


ndim=2; neq=MUA.Nnodes;

hnod=reshape(h(MUA.connectivity,1),MUA.Nele,MUA.nod);
txzbnod=reshape(txzb(MUA.connectivity,1),MUA.Nele,MUA.nod);
tyzbnod=reshape(tyzb(MUA.connectivity,1),MUA.Nele,MUA.nod);


if ~CtrlVar.AGlenisElementBased 
    AGlennod=reshape(AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
    nnod=reshape(n(MUA.connectivity,1),MUA.Nele,MUA.nod);
end

[points,weights]=sample('triangle',MUA.nip,ndim);


bx=zeros(MUA.Nele,MUA.nod);
by=zeros(MUA.Nele,MUA.nod);

MLC=BCs2MLC(CtrlVar,MUA,BCs); Ludvd=MLC.udvdL ; Ludvdrhs=MLC.udvdRhs;

% vector over all elements for each integration point
for Iint=1:MUA.nip
    
    fun=shape_fun(Iint,ndim,MUA.nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    
    if isfield(MUA,'Deriv') && isfield(MUA,'DetJ')
        %Deriv=MUA.Deriv(:,:,:,Iint);
        detJ=MUA.DetJ(:,Iint);
    else
        [~,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
    end
    % Deriv : MUA.Nele x dof x nod
    %  detJ : MUA.Nele
    
    % values at integration point
    
    
    hint=hnod*fun;
    txzbint=txzbnod*fun;
    tyzbint=tyzbnod*fun;
    
    if ~CtrlVar.AGlenisElementBased 
        AGlen=AGlennod*fun;
        n=nnod*fun;
    end
 
    detJw=detJ*weights(Iint);
    D=(2./(n+1)).*AGlen.*(txzbint.^2+tyzbint.^2).^((n-1)/2).*hint.*detJw;  % all variables defined at this integration point
    
    for Inod=1:MUA.nod

        rhs=D.*fun(Inod);
        bx(:,Inod)=bx(:,Inod)+rhs.*txzbint;
        by(:,Inod)=by(:,Inod)+rhs.*tyzbint;
        
    end
end

% assemble right-hand side

rhx=sparseUA(neq,1); rhy=sparseUA(neq,1);
for Inod=1:MUA.nod
    rhx=rhx+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),bx(:,Inod),neq,1);
    rhy=rhy+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),by(:,Inod),neq,1);
end

M=MassMatrix2D1dof(MUA);

if isempty(Ludvd)
    sol=M\[rhx rhy] ;
    ud=full(sol(:,1)) ; vd=full(sol(:,2));
else

    Z=sparse(MUA.Nnodes,MUA.Nnodes);
    M=[ M Z ; Z M ];
    udvdLambda=zeros(size(Ludvd,1),1) ;
    [sol,udvdLambda]=solveKApeSymmetric(M,Ludvd,[rhx ; rhy],Ludvdrhs,[],udvdLambda,CtrlVar);
    ud=full(sol(1:MUA.Nnodes));
    vd=full(sol(MUA.Nnodes+1:2*MUA.Nnodes));

    
end

if ~isreal(ud) ; save TestSave ; error('uvSSHEETplus:udNotReal','ud not real!') ; end 
if ~isreal(vd) ; save TestSave ; error('uvSSHEETplus:vdNotReal','vb not real!') ; end 


end