function [ud,vd]=uvSSHEET(CtrlVar,MUA,BCs,AGlen,n,rho,g,s,h)


%  calculates deformational velocity based on the SSHEET (SIA) approximation
%
%  u and v are nodal velocities
%
%  u= -(2A/(n+1)) (rho g)^n | grad_{xy} s|^(n-1) h^(n+1) \p_x s
%
%  N_p N_q u_q = -N_p (2A/(n+1)) (rho g)^n | grad_{xy} s|^(n-1) h^(n+1) \p_x s
%

narginchk(9,9)

ndim=2; neq=MUA.Nnodes;

hnod=reshape(h(MUA.connectivity,1),MUA.Nele,MUA.nod);
snod=reshape(s(MUA.connectivity,1),MUA.Nele,MUA.nod);
rhonod=reshape(rho(MUA.connectivity,1),MUA.Nele,MUA.nod);

if ~CtrlVar.AGlenisElementBased
    AGlennod=reshape(AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
    nnod=reshape(n(MUA.connectivity,1),MUA.Nele,MUA.nod);
end

[points,weights]=sample('triangle',MUA.nip,ndim);
bx=zeros(MUA.Nele,MUA.nod);
by=zeros(MUA.Nele,MUA.nod);

MLC=BCs2MLC(MUA,BCs);
Ludvd=MLC.udvdL ; Ludvdrhs=MLC.udvdRhs;



% vector over all elements for each integartion point
for Iint=1:MUA.nip
    
    fun=shape_fun(Iint,ndim,MUA.nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    
    if isfield(MUA,'Deriv') && isfield(MUA,'DetJ')
        Deriv=MUA.Deriv(:,:,:,Iint);
        detJ=MUA.DetJ(:,Iint);
    else
        [Deriv,detJ]=derivVector(MUA.coordinates,MUA.connectivity,MUA.nip,Iint);
    end
    % Deriv : MUA.Nele x dof x nod
    %  detJ : MUA.Nele
    
    % values at integration point
    
    
    hint=hnod*fun;
    rhoint=rhonod*fun;
    
    if ~CtrlVar.AGlenisElementBased
        AGlen=AGlennod*fun;
        n=nnod*fun;
    end
    
    
    
    dsdx=zeros(MUA.Nele,1); dsdy=zeros(MUA.Nele,1);
    
    % derivatives for all elements at this integration point
    for Inod=1:MUA.nod
        dsdx=dsdx+Deriv(:,1,Inod).*snod(:,Inod);
        dsdy=dsdy+Deriv(:,2,Inod).*snod(:,Inod);
    end
    
    gradSurf=sqrt(abs(dsdx.*dsdx+dsdy.*dsdy));
    detJw=detJ*weights(Iint);
    
    %T=gradSurf.^(n-1).*(rhoint.*g.*hint).^n;
    
    D=2*AGlen.*gradSurf.^(n-1).*hint.^(n+1).*(rhoint.*g).^n./(n+1);
    
    for Inod=1:MUA.nod
        
        
        rhs=D.*fun(Inod);
        rhsx=-rhs.*dsdx;
        rhsy=-rhs.*dsdy;
        
        bx(:,Inod)=bx(:,Inod)+rhsx.*detJw;
        by(:,Inod)=by(:,Inod)+rhsy.*detJw;
        
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
    sol=M\[rhx rhy] ;  % solve this for two right-hand sides
    ud=full(sol(:,1)) ; vd=full(sol(:,2));
else
    
    Z=sparse(MUA.Nnodes,MUA.Nnodes);
    M=[ M Z ; Z M ];
    udvdLambda=zeros(size(Ludvd,1),1) ;
    [sol,udvdLambda]=solveKApeSymmetric(M,Ludvd,[rhx ; rhy],Ludvdrhs,[],udvdLambda,CtrlVar);
    ud=full(sol(1:MUA.Nnodes));
    vd=full(sol(MUA.Nnodes+1:2*MUA.Nnodes));
    
end




end