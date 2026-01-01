

function K=dFuvdC(CtrlVar,MUA,F)

%%
%
% assembles the matrix K which is the FE form of 
% 
% $$d_{\mathbf{C}} \mathbf{F}_{\mathbf{v}}$$
%
% where
% 
% $$\mathbf{F}_{\mathbf{v}}$$
%
% is the forward model.
% 
% Ruv=Tint-Fext;
% Tint   : internal nodal forces
% Fint   : external nodal forces
%
%
% \p F/\p_ C = [ \p F_1^x/\p C_1    \p F_1^x/\p C_2   ...  \p F_1^x/\p C_n ]
%              [ \p F_2^x/\p C_1    \p F_2^x/\p C_2   ...  \p F_1^x/\p C_n ]
%              [         .                 .                      .        ]
%              [ \p F_n^x/\p C_1    \p F_n^x/\p C_2   ...  \p F_n^x/\p C_n ]
%              [ \p F_1^y/\p C_1    \p F_1^y/\p C_2   ...  \p F_1^y/\p C_n ]
%              [ \p F_2^y/\p C_1    \p F_2^y/\p C_2   ...  \p F_2^y/\p C_n ]
%              [         .                 .                      .        ]
%              [ \p F_n^y/\p C_1    \p F_n^y/\p C_2   ...  \p F_n^y/\p C_n ]
%
%
% $$ F_x=\partial_x ( h \eta ( 4 \partial_x u + 2 \partial_y v)) + \partial_y ( h \eta (\partial_y u + \partial_x v) ) - t_x -   \frac{1}{2} g \partial_x (\rho h^2 -  \rho_o d^2)- g\,\mathcal{H}(h-h_f) (\rho h -\rho_o H^{+}) \partial_x B =0 $$
%
% $$ F_y = \partial_y ( h \eta ( 4 \partial_y v + 2 \partial_x u)) + \partial_x ( h \eta (\partial_x v + \partial_y y) ) - t_y -   \frac{1}{2} g \partial_y (\rho h^2 -  \rho_o d^2)- g\,\mathcal{H}(h-h_f) (\rho h -\rho_o H^{+}) \partial_y B =0 $$
% 
% with
%
% $$ t_x = t_x(C,u,v) $$
%
% $$ t_y = t_y(C,u,v) $$
%
% and
%
% $$ \eta=\eta(A,u,v) $$
%
%
% $$ \partial F_x /\partial C = \partial t_x/\partial C $$
%
%
% The function BasalDrag.m can return a quantity dFuvdC from which we can calculate the partial derivatives as 
%
% dF_x/dC = dFuvdC  u
% dF_y/dC = dFuvdC  v
%

%
%%


ndim=2; dof=2;
nNodes=MUA.Nnodes ;

hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);
Cnod=reshape(F.C(MUA.connectivity,1),MUA.Nele,MUA.nod);
mnod=reshape(F.m(MUA.connectivity,1),MUA.Nele,MUA.nod);

if CtrlVar.SlidingLaw=="Budd"
    qnod=reshape(F.q(MUA.connectivity,1),MUA.Nele,MUA.nod);
else
    qnod=mnod*0 ;  % just to avoid asking this again within a loop
end

if ~isempty(F.muk)
    muknod=reshape(F.muk(MUA.connectivity,1),MUA.Nele,MUA.nod);
else
    muknod=mnod*0 ;
end

if ~isempty(F.V0)
    V0nod=reshape(F.V0(MUA.connectivity,1),MUA.Nele,MUA.nod);
else
    V0nod=mnod*0 ;
end


Bnod=reshape(F.B(MUA.connectivity,1),MUA.Nele,MUA.nod);
Snod=reshape(F.S(MUA.connectivity,1),MUA.Nele,MUA.nod);
rhonod=reshape(F.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);


dFxdC=zeros(MUA.Nele,MUA.nod,MUA.nod);
dFydC=zeros(MUA.Nele,MUA.nod,MUA.nod);

for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    detJ=MUA.DetJ(:,Iint);

    hint=hnod*fun;
    uint=unod*fun;
    vint=vnod*fun;
    Cint=Cnod*fun; Cint(Cint<CtrlVar.Cmin)=CtrlVar.Cmin;
    mint=mnod*fun;
    qint=qnod*fun;
    mukint=muknod*fun;
    V0int=V0nod*fun;
    Bint=Bnod*fun;
    Sint=Snod*fun;
    Hint=Sint-Bint;
    rhoint=rhonod*fun;

    hfint=F.rhow*Hint./rhoint;                                   % this is linear, so fine to evaluate at int in this manner
    Heint = HeavisideApprox(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);  % important to calculate Heint and deltaint in a consistent manner





    %
    % dF/dC=dtaux/dC uAdjoint + dtauy/dC vAdjoint
    %
    % dtaux/dC= He u * dbeta2/dC
    %
    % beta2= (C+CtrlVar.Czero).^(-1./m).*(sqrt(ub.*ub+vb.*vb+CtrlVar.SpeedZero^2)).^(1./m-1) ;
    %

    % setting this CtrlVar field to true ensures that BasalDrag.m returns the (point) derivative
    CtrlVar.Inverse.dFuvdClambda=true;
    Ctemp= ...
        BasalDrag(CtrlVar,MUA,Heint,[],hint,Bint,Hint,rhoint,F.rhow,uint,vint,Cint,mint,[],[],[],[],[],[],[],[],qint,F.g,mukint,V0int);
    CtrlVar.Inverse.dFuvdClambda=false;

    detJw=detJ*MUA.weights(Iint);
    for Inod=1:MUA.nod
        for Jnod=1:MUA.nod

            % \p F_x/C = Ctemp * u
            % \p F_y/C = Ctemp * y

            %  Fx_i =taux.*fun(Inod);
            % \p Fx_i/\p C_j = (\p taux/\p C_j) fun(Inod)  = Ctemp uint fun(Inod) fun(Jnod)

            dFxdC(:,Inod,Jnod)=dFxdC(:,Inod,Jnod)+Ctemp.*uint.*fun(Inod).*fun(Jnod).*detJw;
            dFydC(:,Inod,Jnod)=dFydC(:,Inod,Jnod)+Ctemp.*vint.*fun(Inod).*fun(Jnod).*detJw;

        end
    end
end


Iind=zeros(MUA.nod*MUA.nod*MUA.Nele*2,1,'uint32'); 
Jind=zeros(MUA.nod*MUA.nod*MUA.Nele*2,1,'uint32');

Xval=zeros(MUA.nod*MUA.nod*MUA.Nele*2,1);
istak=0;

for Inod=1:MUA.nod
    %istak=0;

    for Jnod=1:MUA.nod

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=dFxdC(:,Inod,Jnod);
        istak=istak+MUA.Nele;

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+nNodes; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=dFydC(:,Inod,Jnod);
        istak=istak+MUA.Nele;
    

    end
    %K=K+sparse(Iind,Jind,Xval,neq,neq);
end

% tSparse=tic;

K=sparseUA(Iind,Jind,Xval,2*nNodes,nNodes);
% tSparse=toc(tSparse);










end







