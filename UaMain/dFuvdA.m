function K=dFuvdA(CtrlVar,MUA,F)

%%
%
% assembles the matrix K which is the FE form of 
% 
% $$d_{\mathbf{A}} \mathbf{F}_{\mathbf{v}}$$
%
% where
% 
% $$\mathbf{F}_{\mathbf{v}}$$
%
% is the forward model.
% 
%
%
% \p F/\p_A = [ \p F_1^x/\p A_1    \p F_1^x/\p A_2   ...  \p F_1^x/\p A_n ]
%              [ \p F_2^x/\p A_1    \p F_2^x/\p A_2   ...  \p F_1^x/\p A_n ]
%              [         .                 .                      .        ]
%              [ \p F_n^x/\p A_1    \p F_n^x/\p A_2   ...  \p F_n^x/\p A_n ]
%              [ \p F_1^y/\p A_1    \p F_1^y/\p A_2   ...  \p F_1^y/\p A_n ]
%              [ \p F_2^y/\p A_1    \p F_2^y/\p A_2   ...  \p F_2^y/\p A_n ]
%              [         .                 .                      .        ]
%              [ \p F_n^y/\p A_1    \p F_n^y/\p A_2   ...  \p F_n^y/\p A_n ]
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
% $$ \partial F^x_i/\partial A_j =  \langle h (\partial \eta/\partial A)  ( 4 \partial_x u + 2 \partial_y v))\, \phi_j | \partial_x \phi_i \rangle 
% + \langle h (\partial \eta/\partial A) (\partial_y u + \partial_x v)\, \phi_j | \partial_y \phi_i \rangle $$
%
% $$ \partial F^y_i/\partial A_j =  \langle h (\partial \eta/\partial A)  ( 4 \partial_y v + 2 \partial_x u))\, \phi_j | \partial_y \phi_i \rangle 
% + \langle h (\partial \eta/\partial A) (\partial_x v + \partial_y u)\, \phi_j | \partial_x \phi_i \rangle $$
%%


ndim=2; 
nNodes=MUA.Nnodes ;

hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);   % Nele x nod
ubnod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vbnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

Anod=reshape(F.AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
nnod=reshape(F.n(MUA.connectivity,1),MUA.Nele,MUA.nod);

dFxdA=zeros(MUA.Nele,MUA.nod,MUA.nod);
dFydA=zeros(MUA.Nele,MUA.nod,MUA.nod);

for Iint=1:MUA.nip

    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ;
    Deriv=MUA.Deriv(:,:,:,Iint);  % Deriv at integration points
    detJ=MUA.DetJ(:,Iint);

    hint=hnod*fun;

    Aint=Anod*fun;
    Aint(Aint<CtrlVar.AGlenmin)=CtrlVar.AGlenmin;
    nint=nnod*fun;
 
    exx=zeros(MUA.Nele,1);
    eyy=zeros(MUA.Nele,1);
    exy=zeros(MUA.Nele,1);

    for Inod=1:MUA.nod

        exx=exx+Deriv(:,1,Inod).*ubnod(:,Inod);
        eyy=eyy+Deriv(:,2,Inod).*vbnod(:,Inod);
        exy=exy+0.5*(Deriv(:,1,Inod).*vbnod(:,Inod) + Deriv(:,2,Inod).*ubnod(:,Inod));

    end


    [~,~,~,dEtadA]=EffectiveViscositySSTREAM(CtrlVar,Aint,nint,exx,eyy,exy) ;


    detJw=detJ*MUA.weights(Iint);
    for Inod=1:MUA.nod
        for Jnod=1:MUA.nod

            % \p F_x/C = Ctemp * u
            % \p F_y/C = Ctemp * y

            %  Fx_i =taux.*fun(Inod);
            % \p Fx_i/\p C_j = (\p taux/\p C_j) fun(Inod)  = Ctemp uint fun(Inod) fun(Jnod)


            dFxdA(:,Inod,Jnod)=dFxdA(:,Inod,Jnod)+hint.*dEtadA.*((4*exx+2*eyy).*Deriv(:,1,Inod)+2.*exy.*Deriv(:,2,Inod)).*fun(Jnod).*detJw;
            dFydA(:,Inod,Jnod)=dFydA(:,Inod,Jnod)+hint.*dEtadA.*((4*eyy+2*exx).*Deriv(:,2,Inod)+2.*exy.*Deriv(:,1,Inod)).*fun(Jnod).*detJw;

       

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

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod); Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=dFxdA(:,Inod,Jnod);
        istak=istak+MUA.Nele;

        Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod)+nNodes; Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod); Xval(istak+1:istak+MUA.Nele)=dFydA(:,Inod,Jnod);
        istak=istak+MUA.Nele;

    end
end

K=sparseUA(Iind,Jind,Xval,2*nNodes,nNodes);

end





