function K=dFuvdB(CtrlVar,MUA,F)

%%
%
% assembles the matrix K which is the FE form of 
% 
% $$d_{\mathbf{B}} \mathbf{F}_{\mathbf{v}}$$
%
% where
% 
% $$\mathbf{F}_{\mathbf{v}}$$
%
% is the forward model.
% 
%
% $$\left[\begin{array}{cccc} 
% \partial F^x_1 /\partial B_1  & \partial F^x_1 /\partial B_2  & \ldots & \partial F^x_1 /\partial B_n  \\  
% \partial F^x_2 /\partial B_1  & \partial F^x_2 /\partial B_2  & \ldots & \partial F^x_2 /\partial B_n  \\  
%              \vdots              &              \vdots              &  \vdots  &    \vdots              \\  
% \partial F^y_1 /\partial B_1  & \partial F^y_1 /\partial B_2 & \ldots & \partial F^y_1 /\partial B_n  \\
% \partial F^y_2 /\partial B_1  & \partial F^y_2 /\partial B_2 & \ldots & \partial F^y_2 /\partial B_n  \\
%              \vdots              &              \vdots              &  \vdots  &    \vdots              \\  
% \end{array}\right] $$
%
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
%
% Currently only done for grounded area where $B=b$, but once this works for grounded area, will do for the more general
% case...
%
% $$b=\mathcal{G} B + (1-\mathcal{G}) (S-\rho h/\rho_0) $$
%
% The above equation for $b$ contains $h$ on the right hand side as well through $h=s-b$. This is OK if we consider $h$ as
% given. Here we consider the upper surface $s$ as given (i.e. as a input) and therefore solve for $b$, arriving at
%
% $$ b=\frac{1}{1-(1-\mathcal{G}) \rho/\rho_o} \, ( \mathcal{G} B + (1-\mathcal{G}) ( S - \rho s/\rho_o ) ) $$
%
% However, this expression still contains a dependency on $b$ on the right-hand side because $\mathcal{G}$ depends on $h$ and
% therefore on $b$.
%
%
% $$h=s-b$$
%
% $$H=S-B$$
%
% $$h_f=\rho_o H/\rho$$
%
% $$\mathcal{G}=\mathcal{H}(h-h_f)$$
%
% $$H^{+}=\mathcal{H}(H) \, H $$
%
% $$
% d=\mathcal{H}(H) \, (S-b) 
% = \mathcal{H}(h_f-h) \rho h /\rho_o + \mathcal{H}(h-h_f) \; \mathcal{H}(H) \, H 
% = \mathcal{H}(h_f-h) \rho h /\rho_o + \mathcal{H}(h-h_f) \; H^{+}
% = (1-\mathcal{G}) \rho h /\rho_o + \mathcal{G} \; H^{+}
% $$
%
%
% Derivatives for $\mathcal{G}=1$ and $\mathrm{d}\mathcal{G}/\mathrm{d}B=0$ with $s$ and $S$ given:
%
% $$\mathrm{d}s/\mathrm{d}B=\mathrm{d}S/\mathrm{d}B=0$$
%
% $$\mathrm{d}b/\mathrm{d}B=1$$
%
% $$\mathrm{d}H/\mathrm{d}B=-1 $$
%
% $$\mathrm{d} h_f/\mathrm{d}B=(\rho_o/\rho) \; \mathrm{d}H/\mathrm{d}B = -\frac{\rho_0}{\rho} $$
%
% $$ \mathrm{d}\mathcal{G}/\mathrm{d}B = 0 $$
%
% $$ \mathrm{d} H^{+} / \mathrm{d}B = \delta(H) \; \mathrm{d}H/\mathrm{d}B \; H + \mathcal{H}(H) \; \mathrm{d}H/\mathrm{d}B = -\delta(H)\, H - \mathcal{H}(H) $$
%
% $$ \mathrm{d} d/\mathrm{d} B = \mathrm{d} H^{+}/\mathrm{d} B = -\delta(H)\, H - \mathcal{H}(H)$$
%
%
%%


       hfint=F.rhow*Hint./rhoint;                                   % this is linear, so fine to evaluate at int in this manner
            Heint = HeavisideApprox(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);  % important to calculate Heint and deltaint in a consistent manner
            HEint = HeavisideApprox(CtrlVar.kH,hfint-hint,CtrlVar.Hh0);

            deltaint=DiracDelta(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);       % i.e. deltaint must be the exact derivative of Heint
            %Deltaint=DiracDelta(CtrlVar.kH,hfint-hint,CtrlVar.Hh0);      %  although delta is an even function...

            Hposint = HeavisideApprox(CtrlVar.kH,Hint,CtrlVar.Hh0).*Hint;

            % Here we apply the definition of d directly at integration points
            dint=HEint.*rhoint.*hint/F.rhow + Heint.*Hposint ;  % definition of d, applied directly at integration points

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

if isempty(MUA.Deriv)
    [MUA.Deriv,MUA.DetJ]=CalcMuaMeshDerivatives(CtrlVar,MUA);
end


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

            dFxdA(:,Inod,Jnod)=dFxdA(:,Inod,Jnod)+hint.*dEtadA.*((4*exx+2*eyy).*Deriv(:,1,Inod)+2.*exy.*Deriv(:,2,Inod)).*fun(Jnod).*detJw;
            dFydA(:,Inod,Jnod)=dFydA(:,Inod,Jnod)+hint.*dEtadA.*((4*eyy+2*exx).*Deriv(:,2,Inod)+2.*exy.*Deriv(:,1,Inod)).*fun(Jnod).*detJw;





            t1=-F.g*    (rhoint.*hint-F.rhow*dint).*dbdx.*fun(Inod)*ca+ rhoint.*F.g.*hint.*sa.*fun(Inod);
            t2=0.5*F.g.*ca*(rhoint.*hint.^2-F.rhow.*dint.^2).*Deriv(:,1,Inod);
            t3=hint.*etaint.*(4*exx+2*eyy).*Deriv(:,1,Inod);
            t4=hint.*etaint.*2.*exy.*Deriv(:,2,Inod);
            t5=taux.*fun(Inod);

            Tx(:,Inod)=Tx(:,Inod)+(t3+t4+t5).*detJw;
            Fx(:,Inod)=Fx(:,Inod)+(t1+t2).*detJw;



            t1=-F.g*(rhoint.*hint-F.rhow*dint).*dbdy.*fun(Inod)*ca;
            t2=0.5*ca*F.g.*(rhoint.*hint.^2-F.rhow.*dint.^2).*Deriv(:,2,Inod);
            t3=hint.*etaint.*(4*eyy+2*exx).*Deriv(:,2,Inod);
            t4=hint.*etaint.*2.*exy.*Deriv(:,1,Inod);
            t5=tauy.*fun(Inod);

            Ty(:,Inod)=Ty(:,Inod)+(t3+t4+t5).*detJw;
            Fy(:,Inod)=Fy(:,Inod)+(t1+t2).*detJw;













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





