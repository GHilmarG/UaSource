


function dFdhlambda=dIdbq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y,dhdp,dbdp,dBdp)
        

%% Calculates the vector quantity:
%
%
% $$ \langle  \delta_{b_i} F^x \phi_i | \Psi_x \rangle + \langle  \delta_{b_i} F^y \phi_i| \Psi_y \rangle $$
%
%
% Note: Here we only consider grounded ice where $B=b$$. This function does not work correctly for sections with floating
% ice. It is assumed that ALL ice is grounded thorough the domain. 
% 
%% Relates to the solution for the velocity components u and v of the system:
%
% 
% 
% $$
% F^x_i= \left \langle  h \eta \, ( 4 \partial_x u + 2 \partial_y v) \vert  \, \partial_x \phi_i \right \rangle
% + \langle   h \eta \, (\partial_y u + \partial_x v)  \vert  \partial_y \phi_i \rangle
% + \langle \mathcal{G} \beta^2\, u , \phi_i \rangle 
%  - \left \langle \frac{1}{2} g \cos(\alpha) \,  (\rho h^2 -  \rho_o d^2)  \Big\vert \partial_x \phi_i \right \rangle
% + \langle g\, \mathcal{G} \, (\rho h -\rho_o H^{+}) \, \partial_x B \vert  \phi_i \rangle  - \langle \rho g \sin(\alpha) \, h  | \phi_i \rangle   =0 
% $$
%
% $$
% F^y_i= \langle  h \eta \, ( 4 \partial_y v + 2 \partial_x u) \vert \partial_y \phi_i \rangle
% +\langle   h \eta \, (\partial_x v + \partial_y u)  \vert \, \partial_x \phi_i \rangle 
% + \langle \mathcal{G} \, \beta^2 \, v \vert  \phi_i \rangle 
%   - \left \langle \frac{1}{2} g \cos(\alpha) \, (\rho h^2 -  \rho_o d^2) \Big|   \, \partial_y \phi_i \right \rangle
% +  \langle g\, \mathcal{G} \, (\rho h -\rho_o H^{+}) \, \partial_y B \vert \phi_i \rangle=0
% $$
%
%
%
% Here we use
% 
% $$g\, \mathcal{G} \,  (\rho h -\rho_o H^{+}) \, \partial_y B =g\, \mathcal{G} \,  (\rho h -\rho_o H^{+}) \, \partial_y b $$
%
% $\mathcal{G}$ is the floating mask, 1 if grounded, 0 if afloat.
%
% $$\mathcal{G}=\mathcal{H}(h-h_f) $$
%
% where $\mathcal{H}$ is the Heaviside step function and
%
%
% $$h_f=(S-B) \rho_o/\rho $$
%
% where:
% 
% $h=s-b$ is the ice thickness
%
% $\rho$ the ice density
%
% $\rho_o$ the ocean density 
%
% $B$ the bedrock
%
% $s$ the upper glacier surface
%
% $b$ the lower glacier surface
%
% $S$ the ocean surface
%
% $$\alpha$$ the slope of the vertical axis of the coordinate system with respect to gravity
%
% $u$ the $x$ velocity component
%
% $v$ the $y$ velocity component
%
% $d$ is the submarine ice thickness (always positive), defined as: 
% 
% $$d=\mathcal{H}(h_f-h) \, \rho h / \rho_o + \mathcal{H}(h-h_f) \, H^{+} $$
%
% which we can also write as
%
% $$d= (1-\mathcal{G} )\, \frac{\rho}{\rho_o}  h + \mathcal{G} \, H^{+} $$
%
% $$H^{+} = \mathcal{H}(H) \, H $$
%
% $$H=S-B $$
%
% The effective viscosity is: 
%
% $$
% \eta= \frac{1}{2} A^{-1/n} \, \left ((\partial_x u)^2 + (\partial_y v)^2 + \partial_x u \,\partial_y v + (\partial_x v + \partial_y u)^2/4+\epsilon_0^2 \right)^{(1-n)/2n} +\eta_0
% $$
%
%
% The effective viscosity is therefore a function of the velocity components and the rheological parameters $A$ and $n$.
%
% The function
% 
%   EffectiveViscositySSTREAM.m
%
% returns the effective viscosity, eta, as well as some derivatives with respect to $A$.
%
% In the particular case of Weertman sliding law $$\beta^2$$ is given by:
%
% $$
% \beta^2=(C+C_0)^{-1/m} \; \left (u_b^2+v_b^2+u_0^2 \right)^{(1-m)/2m} 
% $$
%
% $$\beta^2$$ is therefore a function of the velocity components, and the basal sliding law parameter $C$. For more general sliding
% laws $\beta^2$ may depend on some other parameters as well.  
%
% Often the basal drag term is written as
%
%
% $$t_{bx} =\mathcal{G} \beta^2\, u $$
%
% $$t_{by} =\mathcal{G} \beta^2\, v $$
%
%
% where $t_{bx}$ and $t_{by}$ are the basal traction components.
%
% The function
% 
%   BasalDrag.m
%
% returns $t_{bx}$ and $t_{by}$ as well as various derivatives with respect to $u$, $v$,  $h$ and $C$
%
%%





narginchk(10,10)

ndim=2;


if ~CtrlVar.IncludeMelangeModelPhysics
    uoint=[];
    voint=[];
    Coint=[];
    moint=[];
    uaint=[];
    vaint=[];
    Caint=[];
    maint=[];
else
    error('Inversion with MelangeModelPhysics not implemented.\n')
end

ca=cos(F.alpha); sa=sin(F.alpha);


hnod=reshape(F.h(MUA.connectivity,1),MUA.Nele,MUA.nod);
snod=reshape(F.s(MUA.connectivity,1),MUA.Nele,MUA.nod);
bnod=reshape(F.b(MUA.connectivity,1),MUA.Nele,MUA.nod);
Bnod=reshape(F.B(MUA.connectivity,1),MUA.Nele,MUA.nod);
Snod=reshape(F.S(MUA.connectivity,1),MUA.Nele,MUA.nod);

dbdpnod=reshape(dbdp(MUA.connectivity,1),MUA.Nele,MUA.nod);
dhdpnod=reshape(dhdp(MUA.connectivity,1),MUA.Nele,MUA.nod);
dBdpnod=reshape(dBdp(MUA.connectivity,1),MUA.Nele,MUA.nod);


rhonod=reshape(F.rho(MUA.connectivity,1),MUA.Nele,MUA.nod);

unod=reshape(F.ub(MUA.connectivity,1),MUA.Nele,MUA.nod);
vnod=reshape(F.vb(MUA.connectivity,1),MUA.Nele,MUA.nod);

Psi_x_node=reshape(Psi_x(MUA.connectivity,1),MUA.Nele,MUA.nod);
Psi_y_node=reshape(Psi_y(MUA.connectivity,1),MUA.Nele,MUA.nod);

AGlennod=reshape(F.AGlen(MUA.connectivity,1),MUA.Nele,MUA.nod);
nnod=reshape(F.n(MUA.connectivity,1),MUA.Nele,MUA.nod);

Cnod=reshape(F.C(MUA.connectivity,1),MUA.Nele,MUA.nod);
mnod=reshape(F.m(MUA.connectivity,1),MUA.Nele,MUA.nod);


if ~isempty(F.q)
    qnod=reshape(F.q(MUA.connectivity,1),MUA.Nele,MUA.nod);
end

if ~isempty(F.V0)
    V0nod=reshape(F.V0(MUA.connectivity,1),MUA.Nele,MUA.nod);
end


if ~isempty(F.muk)
    muknod=reshape(F.muk(MUA.connectivity,1),MUA.Nele,MUA.nod);
end


% [points,weights]=sample('triangle',MUA.nip,ndim);
T=zeros(MUA.Nele,MUA.nod);


hfnod=F.rhow*(Snod-Bnod)./rhonod;

for Iint=1:MUA.nip
    
    fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    detJ=MUA.DetJ(:,Iint);
    Deriv=MUA.Deriv(:,:,:,Iint);
    
    
    hint=hnod*fun;
    bint=bnod*fun;
    
    uint=unod*fun;
    vint=vnod*fun;

    rhoint=rhonod*fun;
    nint=nnod*fun;
    Cint=Cnod*fun; Cint(Cint<CtrlVar.Cmin)=CtrlVar.Cmin;
    mint=mnod*fun;
    Bint=Bnod*fun;
    Sint=Snod*fun;
    Hint=Sint-Bint;


    AGlenInt=AGlennod*fun;
    AGlenInt(AGlenInt<CtrlVar.AGlenmin)=CtrlVar.AGlenmin;

    if ~isempty(F.q)
        qint=qnod*fun;
    else
        qint=[];
    end

    if ~isempty(F.muk)
        mukint=muknod*fun;
    else
        mukint=[];
    end

    if ~isempty(F.V0)
        V0int=V0nod*fun;
    else
        V0int=[];
    end


    Psi_x_int=Psi_x_node*fun;
    Psi_y_int=Psi_y_node*fun;
    
    hfint=F.rhow*Hint./rhoint;
    
    % hfint=hfnod*fun;
    
    deltaint=DiracDelta(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);
    Heint = HeavisideApprox(CtrlVar.kH,hint-hfint,CtrlVar.Hh0);
    
    HeHint = HeavisideApprox(CtrlVar.kH,Hint,CtrlVar.Hh0);
    deltaHint=DiracDelta(CtrlVar.kH,Hint,CtrlVar.Hh0);
    dint = HeHint.*(Sint-bint);  % draft
    %dint = HeavisideApprox(CtrlVar.kH,Hint,CtrlVar.Hh0).*(Sint-bint);  % draft
    
    dhdpint=dhdpnod*fun;
    dbdpint=dbdpnod*fun;
    dBdpint=dBdpnod*fun;
    
    % only correct for B !!  +gera+
    dhdpint= -Heint ;
    
    dbdpint= Heint ;
    dBdpint=1 ;
    
    
    dlxdx=zeros(MUA.Nele,1);    dlydx=zeros(MUA.Nele,1);    dlxdy=zeros(MUA.Nele,1);    dlydy=zeros(MUA.Nele,1);
    dsdx=zeros(MUA.Nele,1);     dsdy=zeros(MUA.Nele,1);
    %dbdx=zeros(MUA.Nele,1);     dbdy=zeros(MUA.Nele,1);
    dhdx=zeros(MUA.Nele,1);     dhdy=zeros(MUA.Nele,1);
    dhfdx=zeros(MUA.Nele,1);     dhfdy=zeros(MUA.Nele,1);
    
    exx=zeros(MUA.Nele,1);
    eyy=zeros(MUA.Nele,1);
    exy=zeros(MUA.Nele,1);
    
    ddbdpdx=zeros(MUA.Nele,1);
    ddbdpdy=zeros(MUA.Nele,1);
    ddhdpdx=zeros(MUA.Nele,1);
    ddhdpdy=zeros(MUA.Nele,1);
    
    for Inod=1:MUA.nod
        
        dsdx=dsdx+Deriv(:,1,Inod).*snod(:,Inod);
        dsdy=dsdy+Deriv(:,2,Inod).*snod(:,Inod);
        
        dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
        dhdy=dhdy+Deriv(:,2,Inod).*hnod(:,Inod);
        
        dhfdx=dhfdx+Deriv(:,1,Inod).*hfnod(:,Inod);
        dhfdy=dhfdy+Deriv(:,2,Inod).*hfnod(:,Inod);
        
        exx=exx+Deriv(:,1,Inod).*unod(:,Inod);
        eyy=eyy+Deriv(:,2,Inod).*vnod(:,Inod);
        exy=exy+0.5*(Deriv(:,1,Inod).*vnod(:,Inod) + Deriv(:,2,Inod).*unod(:,Inod));
        
        
        dlxdx=dlxdx+Deriv(:,1,Inod).*Psi_x_node(:,Inod);
        dlxdy=dlxdy+Deriv(:,2,Inod).*Psi_x_node(:,Inod);
        
        dlydx=dlydx+Deriv(:,1,Inod).*Psi_y_node(:,Inod);
        dlydy=dlydy+Deriv(:,2,Inod).*Psi_y_node(:,Inod);
        
        ddbdpdx=ddbdpdx+Deriv(:,1,Inod).*dbdpnod(:,Inod);
        ddbdpdy=ddbdpdy+Deriv(:,2,Inod).*dbdpnod(:,Inod);
        
        ddhdpdx=ddhdpdx+Deriv(:,1,Inod).*dhdpnod(:,Inod);
        ddhdpdy=ddhdpdy+Deriv(:,2,Inod).*dhdpnod(:,Inod);
        
    end
    
    dbdx=dsdx-dhdx; dbdy=dsdy-dhdy;
    
    detJw=detJ*MUA.weights(Iint);
    
    
    
    
    %
    %  dh/db= d(s-b)/db = ds/db - db/db = ds/db -1 =
    %
    % Using: s= (1-rhow/rho) b + rhow S/rho   for h<h_f
    %
    % I assume that s is independent of b where grounded, i.e. ds/db = 0,
    % on the other hand I can not assume that s is independent of b where afloat because
    % that will violate the floating condition.
    %
    % Therefore:  ds/db =  (1-Heint) (1-rhow/rho)
    %     %
    %     % dh/db = (1-Heint) (1-rhow/rho) -1
    %     if contains(CtrlVar.Inverse.InvertFor,'-B-')
    %         % only change B, i.e. dhdb=dhdB*Heing
    %         dhdbint= -Heint;
    %     else
    %
    %         dhdbint= (1-Heint).*(1-F.rhow./rhoint)-1 ;
    %     end
    %     %    dhdb=-1;
    
    
  
    [~,~,~,~,~,~,dtaubxdh,dtaubydh] = BasalDrag(CtrlVar,[],Heint,deltaint,hint,Bint,Hint,rhoint,F.rhow,uint,vint,Cint,mint,...
        uoint,voint,Coint,moint,uaint,vaint,Caint,maint,...
        qint,F.g,mukint,V0int);
    etaint=EffectiveViscositySSTREAM(CtrlVar,AGlenInt,nint,exx,eyy,exy);
    
    
    for Inod=1:MUA.nod
        
        
        % uvMatrixAssembly:
        %
        %         t1=-F.g*(rhoint.*hint-F.rhow*dint).*dbdx.*fun(Inod)*ca+ rhoint.*F.g.*hint.*sa.*fun(Inod);
        %         t2=0.5*F.g.*ca*(rhoint.*hint.^2-F.rhow.*dint.^2).*Deriv(:,1,Inod);
        %
        %         t3=hint.*etaint.*(4*exx+2*eyy).*Deriv(:,1,Inod);
        %         t4=hint.*etaint.*2.*exy.*Deriv(:,2,Inod);
        %         t5=taux.*fun(Inod); % beta2int.*uint.*fun(Inod);  % basal friction, Weertman, u
        %
        %         Tx(:,Inod)=Tx(:,Inod)+(t3+t4+t5).*detJw;
        %         Fx(:,Inod)=Fx(:,Inod)+(t1+t2).*detJw;
        %
        %         t1=-F.g*ca*(rhoint.*hint-F.rhow*dint).*dbdy.*fun(Inod);
        %         t2=0.5*ca*F.g.*(rhoint.*hint.^2-F.rhow.*dint.^2).*Deriv(:,2,Inod);
        %
        %         t3=hint.*etaint.*(4*eyy+2*exx).*Deriv(:,2,Inod);
        %         t4=hint.*etaint.*2.*exy.*Deriv(:,1,Inod);
        %         t5=tauy.*fun(Inod);                       % beta2int.*vint.*fun(Inod); % basal friction, Weertman, v
        %
        %         Ty(:,Inod)=Ty(:,Inod)+(t3+t4+t5).*detJw;
        %         Fy(:,Inod)=Fy(:,Inod)+(t1+t2).*detJw;
        %
        %
        
        %dtaubxdh=0;
        %dtaubydh=0;

        %          dbdpint=Heint;
        %          ddbdpdx=deltaint.*(dhdx-dhfdx);
        %          ddbdpdy=deltaint.*(dhdy-dhfdy);
        
        
        % Note: db/dx  needs to be perturbed with respect to B (i.e. p)
        %  b = G B  + (1-G) ... (not function of B)
        % d (delta b)/dx = d (delta (G B) ) / dx
        %                = d (deltaG  B) ) / dx + d (G delta B) / dx
        %                =     ?                + dG/dx delta B + G d(delta B) /dx
        %                =     ?                + dG/dx phi + G d(phi)/dx
        %                =     ?                + dG/dx fun + G deriv
        %                =     ?                + delta(h-h_f) d(h-h_f)/dx     fun + G deriv
        %
        % if db/dB = db/dp = G then  set G=dbdp
        %
        %   t1=-ca*F.g*(rhoint.*hint-F.rhow*dint).*dbdx.*fun(Inod)+ rhoint.*F.g.*hint.*sa.*fun(Inod);
        %t1=   (-ca*F.g* (rhoint.*hint-F.rhow*dint)   .*(ddbdpdx.*fun(Inod)+dbdpint.*Deriv(:,1,Inod)) ...  %der(1,Inod)
        %       -ca*F.g*(rhoint.*dhdpint+F.rhow*HeHint.*dBdpint).*dbdx  .*fun(Inod) ...
        %    + rhoint.*F.g.*sa.*dhdpint.*fun(Inod)                                ).*uAdjointint;
        
        %         t1=-F.g*(rhoint.*hint-F.rhow*dint).*dbdx.*fun(Inod)*ca+ rhoint.*F.g.*hint.*sa.*fun(Inod);
        test1=1;
        test2=1; 
        
        t1=-ca*F.g*...
            (...
              (rhoint.*hint-F.rhow*dint).*(test1*deltaint.*(dhdx-dhfdx).*fun(Inod)+ dbdpint.*Deriv(:,1,Inod))...
              +(rhoint.*dhdpint.*fun(Inod)+F.rhow*(HeHint.*dbdpint+test2*deltaHint.*dBdpint.*(Sint-bint)).*fun(Inod)).*dbdx...
              ).*Psi_x_int ...
            +rhoint.*F.g.*sa.*dhdpint.*fun(Inod).*Psi_x_int;
        
        %         t2=0.5*F.g.*ca*(rhoint.*hint.^2-F.rhow.*dint.^2).*Deriv(:,1,Inod);
        
        % t2=0.5*F.g.*ca*(rhoint.*hint.^2-F.rhow.*dint.^2).*Deriv(:,1,Inod);
        t2=ca*F.g.*(rhoint.*hint.*dhdpint.*fun(Inod)-F.rhow.*dint.*(-HeHint.*dbdpint-test2*deltaHint.*dBdpint.*(Sint-bint)).*fun(Inod)).*dlxdx;
        
        t3=dhdpint.*fun(Inod).*etaint.*(4*exx+2*eyy).*dlxdx;
        t4=dhdpint.*fun(Inod).*etaint.*2.*exy.*dlxdy;
        t5=(dhdpint+F.rhow*dBdpint./rhoint) .*dtaubxdh.*Psi_x_int.*fun(Inod);
        t5=0;
        
        Fx=(t1+t2).*detJw;
        %         t1=-F.g*ca*(rhoint.*hint-F.rhow*dint).*dbdy.*fun(Inod);
        %t1=   (-ca*F.g* (rhoint.*hint-F.rhow*dint)   .*(ddbdpdy.*fun(Inod)+dbdpint.*Deriv(:,2,Inod))...
        %    -ca*F.g*(rhoint.*dhdpint+F.rhow*HeHint.*dBdpint).*dbdy  .*fun(Inod) ...
        %                                                                       ).*vAdjointint;
        
        %         t1=-F.g*ca*(rhoint.*hint-F.rhow*dint).*dbdy.*fun(Inod);
        t1=-F.g*ca*...
            (  (rhoint.*hint-F.rhow*dint).*(test1*deltaint.*(dhdy-dhfdy).*fun(Inod)+dbdpint.*Deriv(:,2,Inod))...
            +(rhoint.*dhdpint.*fun(Inod)+F.rhow*(HeHint.*dbdpint+test2*deltaHint.*dBdpint.*(Sint-bint)).*fun(Inod)).*dbdy)...
            .*Psi_y_int;
        
        Tx=(t3+t4+t5).*detJw;
        
        
        % t2=0.5*ca*g.*(rhoint.*hint.^2-F.rhow.*dint.^2).*Deriv(:,2,Inod);
        t2=F.g*ca*(rhoint.*hint.*dhdpint.*fun(Inod)-F.rhow.*dint.*(-HeHint.*dbdpint-test2*deltaHint.*dBdpint.*(Sint-bint)).*fun(Inod)).*dlydy ; 
        
        t3=dhdpint.*fun(Inod).*etaint.*(4*eyy+2*exx).*dlydy; % t3=hint.*etaint.*(4*eyy+2*exx).*Deriv(:,2,Inod);
        t4=dhdpint.*fun(Inod).*etaint.*2.*exy.*dlydx ; % t4=hint.*etaint.*2.*exy.*Deriv(:,1,Inod);
        t5=(dhdpint+F.rhow*dBdpint./rhoint) .*dtaubydh.*Psi_y_int.*fun(Inod);   % 5=tauy.*fun(Inod);
        t5=0; 
        
        
        
        Fy=(t1+t2).*detJw;
        Ty=(t3+t4+t5).*detJw;
        
        
        T(:,Inod)=T(:,Inod)-Tx+Fx-Ty+Fy;   % opposite sign to K because of the Newton sign
        
    end
    
    
    
end

dFdhlambda=zeros(MUA.Nnodes,1);


for Inod=1:MUA.nod
    dFdhlambda=dFdhlambda+sparse(MUA.connectivity(:,Inod),ones(MUA.Nele,1),T(:,Inod),MUA.Nnodes,1);
end

dFdhlambda=-dFdhlambda;

% dFdhlambda=ApplyAdjointGradientPreMultiplier(CtrlVar,MUA,BCsAdjoint,CtrlVar.Inverse.AdjointGradient.UseBCs.B,dFdhlambda);
% Now this is done for the whole assembled dIdp gradient
 
end






