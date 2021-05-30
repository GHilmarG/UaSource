function [ud,vd,ub,vb]=uvSSHEET(CtrlVar,MUA,BCs,AGlen,n,C,m,rho,g,s,h)


%  calculates deformational and basal velocity based on the SSHEET (SIA) approximation
%
%  u and v are nodal velocities
%
%  u= -(2A/(n+1)) (rho g)^n | grad_{xy} s|^(n-1) h^(n+1) \p_x s
%
%  N_p N_q u_q = -N_p (2A/(n+1)) (rho g)^n | grad_{xy} s|^(n-1) h^(n+1) \p_x s
%

narginchk(11,11)
nargoutchk(4,4)


ndim=2; neq=MUA.Nnodes;

hnod=reshape(h(MUA.connectivity,1),MUA.Nele,MUA.nod);
snod=reshape(s(MUA.connectivity,1),MUA.Nele,MUA.nod);
rhonod=reshape(rho(MUA.connectivity,1),MUA.Nele,MUA.nod);
AGlennod=reshape(AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
nnod=reshape(n(MUA.connectivity,1),MUA.Nele,MUA.nod);
Cnod=reshape(C(MUA.connectivity,1),MUA.Nele,MUA.nod);
mnod=reshape(m(MUA.connectivity,1),MUA.Nele,MUA.nod);

% [points,weights]=sample('triangle',MUA.nip,ndim);

bxd=zeros(MUA.Nele,MUA.nod);
byd=zeros(MUA.Nele,MUA.nod);
bxb=zeros(MUA.Nele,MUA.nod);
byb=zeros(MUA.Nele,MUA.nod);


MLC=BCs2MLC(CtrlVar,MUA,BCs);
Ludvd=MLC.udvdL ; Ludvdrhs=MLC.udvdRhs;
Lubvb=MLC.ubvbL ; Lubvbrhs=MLC.ubvbRhs;



% vector over all elements for each integartion point
for Iint=1:MUA.nip
    
    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    
    
    Deriv=MUA.Deriv(:,:,:,Iint);
    detJ=MUA.DetJ(:,Iint);
    
    % Deriv : MUA.Nele x dof x nod
    %  detJ : MUA.Nele
    
    % values at integration point
    
    
    hint=hnod*fun;
    rhoint=rhonod*fun;
    
    
    AGlen=AGlennod*fun;
    n=nnod*fun;
      
    C=Cnod*fun;
    m=mnod*fun;
    
    
    
    dsdx=zeros(MUA.Nele,1); dsdy=zeros(MUA.Nele,1);
    
    % derivatives for all elements at this integration point
    for Inod=1:MUA.nod
        dsdx=dsdx+Deriv(:,1,Inod).*snod(:,Inod);
        dsdy=dsdy+Deriv(:,2,Inod).*snod(:,Inod);
    end
    
    gradSurf=sqrt(abs(dsdx.*dsdx+dsdy.*dsdy));
    detJw=detJ*MUA.weights(Iint);
    
    %T=gradSurf.^(n-1).*(rhoint.*g.*hint).^n;
    
    Dd=2*AGlen.* (gradSurf.^(n-1)) .* (hint.^(n+1)) .* ((rhoint.*g).^n) ./(n+1);
    
    Db=C.*(gradSurf.^(m-1)) .* ((rhoint.*g.*hint).^m) ;
    
    for Inod=1:MUA.nod
        
        % Deformational
        rhsd=Dd.*fun(Inod);
        rhsxd=-rhsd.*dsdx;
        rhsyd=-rhsd.*dsdy;
        
        bxd(:,Inod)=bxd(:,Inod)+rhsxd.*detJw;
        byd(:,Inod)=byd(:,Inod)+rhsyd.*detJw;
        
        
        % basal
        rhsb=Db.*fun(Inod);
        rhsxb=-rhsb.*dsdx;
        rhsyb=-rhsb.*dsdy;
        
        bxb(:,Inod)=bxb(:,Inod)+rhsxb.*detJw;
        byb(:,Inod)=byb(:,Inod)+rhsyb.*detJw;
        
    end
end

% assemble right-hand side

rhxd=sparseUA(neq,1); rhyd=sparseUA(neq,1);
rhxb=sparseUA(neq,1); rhyb=sparseUA(neq,1);
for Inod=1:MUA.nod
    rhxd=rhxd+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),bxd(:,Inod),neq,1);
    rhyd=rhyd+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),byd(:,Inod),neq,1);
    rhxb=rhxb+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),bxb(:,Inod),neq,1);
    rhyb=rhyb+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),byb(:,Inod),neq,1);
end


%M=MassMatrix2D1dof(MUA);
M=MUA.M; 


% Deformational
if isempty(Ludvd)
    sol=M\[rhxd rhyd] ;  % solve this for two right-hand sides
    ud=full(sol(:,1)) ; vd=full(sol(:,2));
    
    
else
    
    Z=sparse(MUA.Nnodes,MUA.Nnodes);
    MZ=[ M Z ; Z M ];
    udvdLambda=zeros(size(Ludvd,1),1) ;
    sol=solveKApeSymmetric(MZ,Ludvd,[rhxd ; rhyd],Ludvdrhs,[],udvdLambda,CtrlVar);
    ud=full(sol(1:MUA.Nnodes));
    vd=full(sol(MUA.Nnodes+1:2*MUA.Nnodes));
    
end


% Basal
if isempty(Lubvb)
    sol=M\[rhxb rhyb] ;  % solve this for two right-hand sides
    ub=full(sol(:,1)) ; vb=full(sol(:,2));
    
    
else
    
    Z=sparse(MUA.Nnodes,MUA.Nnodes);
    MZ=[ M Z ; Z M ];
    ubvbLambda=zeros(size(Lubvb,1),1) ;
    sol=solveKApeSymmetric(MZ,Lubvb,[rhxb ; rhyb],Lubvbrhs,[],ubvbLambda,CtrlVar);
    ub=full(sol(1:MUA.Nnodes));
    vb=full(sol(MUA.Nnodes+1:2*MUA.Nnodes));
    
end




end