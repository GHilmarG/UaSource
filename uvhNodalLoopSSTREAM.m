function [Tx,Fx,Ty,Fy,Th,Fh,Kxu,Kxv,Kyu,Kyv,Kxh,Kyh,Khu,Khv,Khh]=...
    uvhNodalLoopSSTREAM(detJw,nod,theta,tau0,Ronly,...
    CtrlVar,Tx,Fx,Ty,Fy,Th,Fh,Kxu,Kxv,Kyu,Kyv,Kxh,Kyh,Khu,Khv,Khh, ...
    Deriv,fun,...
    exx,eyy,exy,exx0,eyy0,...
    dhdx,dhdy,dh0dx,dh0dy,drhodx,drhody,dbdx,dbdy,dBdx,dBdy,Hposint,Dddhint,...
    ca,sa,g,dt,...
    etaint,Eint,...
    h0barr,h1barr,...
    taux,tauy,dtauxdu,dtauxdv,dtauydu,dtauydv,dtauxdh,dtauydh,...
    Heint,deltaint,rhoint,rhow,uint,vint,u0int,v0int,dint,...
    hint,h0int,a1int,a0int,dadhint,lambda_h)
% 
% if ~Ronly
%     save uvhNodalLoopSSTREAM.mat
%     
% end

%%
% 
% $$ t_1  = -c_a \, g \,  (\rho \, h - \rho_w d  ) \, db/dx \, N_i  + \rho g \, s_a \, N_i \, h   $$
%
% $$ D_h \{ t_1 \} = -c_a \, g \,  (\rho \, \delta h- \rho_w D_h \{ d \} ) \, db/dx \, N_i  + \rho g \, s_a \, N_i \, \delta h   $$
% 
% $$ d=\mathcal{H}(h_f-h) \, \rho h /\rho_w + \mathcal{H}(h-h_f) \,  H^{+} $$
%
% $$ D_h\{d\} = -\delta(h_f-h) \, \delta h \, \rho h /\rho_w + \mathcal{H}(h_f-h) \, \rho \, \delta h /\rho_w + \delta(h-h_f) \, H^{+} \, \delta h $$
%
% $$ D_h\{d\} = \delta(h_f-h) \, \left ( -\rho h /\rho_w + H^{+})   + \mathcal{H}(h_f-h) \, \rho/\rho_w ) \, \right )\delta h $$
%
% $$ D_h \{ t_1 \} = -c_a \, g \,  (\rho \, \delta h- \rho_w \left ( \delta(h_f-h) \, \left ( -\rho h /\rho_w + H^{+} \right )   + \mathcal{H}(h_f-h) \, \rho/\rho_w ) \, \right  ) \delta h    ) \, db/dx \, N_i  + \rho g \, s_a \, N_i \, \delta h   $$
%
% $$ D_h \{ t_1 \} = - c_a \, g \,  \left ( \rho \, \delta h +  \left (\delta(h_f-h) \, \left ( \rho h - H^{+} \rho_w \right )   - \mathcal{H}(h_f-h) \, \rho \, \right ) \delta h  \right ) \, db/dx \, N_i  + \rho g \, s_a \, N_i \, \delta h   $$
%
% $$ D_h \{ t_1 \} = - c_a \, g \,  \left ( \rho +  \delta(h_f-h) \, \left ( \rho h - H^{+} \rho_w \right )   - \mathcal{H}(h_f-h) \, \rho \, \right ) \, d_xb \, N_i \, \delta h + \rho g \, s_a \, N_i \, \delta h   $$
% 
%%

qx1dx=rhoint.*exx.*hint+rhoint.*uint.*dhdx+drhodx.*uint.*hint;
qy1dy=rhoint.*eyy.*hint+rhoint.*vint.*dhdy+drhody.*vint.*hint;


qx0dx=rhoint.*exx0.*h0int+rhoint.*u0int.*dh0dx+drhodx.*u0int.*h0int;
qy0dy=rhoint.*eyy0.*h0int+rhoint.*v0int.*dh0dy+drhody.*v0int.*h0int;

% Above, corrected from
% qx0dx=rhoint.*exx0.*h0int+rhoint.*u0int.*dh0dx+drhodx.*u0int.*hint;
% qy0dy=rhoint.*eyy0.*h0int+rhoint.*v0int.*dh0dy+drhody.*v0int.*hint;
% on 22/09/2021 and again on 13/03/2022


for Inod=1:nod

    SUPG=fun(Inod)+(1-theta).*tau0.*(u0int.*Deriv(:,1,Inod)+v0int.*Deriv(:,2,Inod));
    funI=fun(Inod) ;
    
    if ~Ronly
        for Jnod=1:nod
            
            Deu=Eint.*((2*exx+eyy).*Deriv(:,1,Jnod)+exy.*Deriv(:,2,Jnod));
            Dev=Eint.*((2*eyy+exx).*Deriv(:,2,Jnod)+exy.*Deriv(:,1,Jnod));
            
            E11=  hint.*(4.*exx+2.*eyy).*Deu.*Deriv(:,1,Inod)+2*hint.*exy.*Deu.*Deriv(:,2,Inod);
            E12=  hint.*(4.*exx+2.*eyy).*Dev.*Deriv(:,1,Inod)+2*hint.*exy.*Dev.*Deriv(:,2,Inod);
            E22=  hint.*(4.*eyy+2.*exx).*Dev.*Deriv(:,2,Inod)+2*hint.*exy.*Dev.*Deriv(:,1,Inod);
            E21=  hint.*(4.*eyy+2.*exx).*Deu.*Deriv(:,2,Inod)+2*hint.*exy.*Deu.*Deriv(:,1,Inod);
            
            
            % Derivative of x-momentum (Fx) with respect to u
            Kxu(:,Inod,Jnod)=Kxu(:,Inod,Jnod)...
                +(4*hint.*etaint.*Deriv(:,1,Inod).*Deriv(:,1,Jnod)...
                +hint.*etaint.*Deriv(:,2,Inod).*Deriv(:,2,Jnod)...
                +E11...
                +dtauxdu.*fun(Jnod).*funI...
                ).*detJw;
            
            
            Kyv(:,Inod,Jnod)=Kyv(:,Inod,Jnod)...
                +(4*hint.*etaint.*Deriv(:,2,Inod).*Deriv(:,2,Jnod)...
                +hint.*etaint.*Deriv(:,1,Inod).*Deriv(:,1,Jnod)...
                +E22...
                +dtauydv.*fun(Jnod).*funI...
                ).*detJw ;
            
            Kxv(:,Inod,Jnod)=Kxv(:,Inod,Jnod)...
                +(etaint.*hint.*(2*Deriv(:,1,Inod).*Deriv(:,2,Jnod)+Deriv(:,2,Inod).*Deriv(:,1,Jnod))...
                +E12...
                +dtauxdv.*fun(Jnod).*funI...
                ).*detJw;


            Kyu(:,Inod,Jnod)=Kyu(:,Inod,Jnod)...
                +(etaint.*hint.*(2*Deriv(:,2,Inod).*Deriv(:,1,Jnod)+Deriv(:,1,Inod).*Deriv(:,2,Jnod))...
                +E21...
                +dtauydu.*fun(Jnod).*funI...    % +Dbeta2Duvint*fun(Jnod).*fun(Inod)...
                ).*detJw;

            % Derivative of Fx momentum with respect to h
            % dFx/dh

            if CtrlVar.Development.Pre2025uvhAssembly


                Kxh(:,Inod,Jnod)=Kxh(:,Inod,Jnod)...
                    +(etaint.*(4*exx+2*eyy).*Deriv(:,1,Inod).*fun(Jnod)...
                    +etaint.*2.*exy.*Deriv(:,2,Inod).*fun(Jnod)...
                    +dtauxdh.*funI.*fun(Jnod)... % +deltaint.*beta2int.*uint.*fun(Inod).*fun(Jnod)..
                    +ca*g*rhoint.*Heint.*dBdx.*funI.*fun(Jnod)...                               % t1:
                    +ca*g*deltaint.*(rhoint.*hint-rhow*Hposint).*dBdx.*funI.*fun(Jnod)... ;     % t1
                    -sa*g*rhoint.*funI.*fun(Jnod)...                                            % t1
                    -ca*g*(rhoint.*hint-rhow*dint.*Dddhint).*Deriv(:,1,Inod).*fun(Jnod)...  ;   % t2
                    ).*detJw;

                Kyh(:,Inod,Jnod)=Kyh(:,Inod,Jnod)...
                    +(etaint.*(4*eyy+2*exx).*Deriv(:,2,Inod).*fun(Jnod)...
                    +etaint.*2.*exy.*Deriv(:,1,Inod).*fun(Jnod)...
                    +dtauydh.*funI.*fun(Jnod)...   % +deltaint.*beta2int.*vint.*fun(Inod).*fun(Jnod)...
                    +ca*g*rhoint.*Heint.*dBdy.*funI.*fun(Jnod)...                                % t1
                    +ca*g*deltaint.*(rhoint.*hint-rhow*Hposint).*dBdy.*funI.*fun(Jnod)... ;      % t1
                    -ca*g*(rhoint.*hint-rhow*dint.*Dddhint).*Deriv(:,2,Inod).*fun(Jnod)...  ;    % t2
                    ).*detJw;


            else


                Kxh(:,Inod,Jnod)=Kxh(:,Inod,Jnod)...
                    +(etaint.*(4*exx+2*eyy).*Deriv(:,1,Inod).*fun(Jnod)...
                    +etaint.*2.*exy.*Deriv(:,2,Inod).*fun(Jnod)...
                    +dtauxdh.*funI.*fun(Jnod)... % +deltaint.*beta2int.*uint.*fun(Inod).*fun(Jnod)..
                    +ca*g*( rhoint - rhow * Dddhint).* dbdx .*funI.*fun(Jnod)...                % t1:
                    -sa*g*rhoint.*funI.*fun(Jnod)...                                            % t1
                    -ca*g*(rhoint.*hint-rhow*dint.*Dddhint).*Deriv(:,1,Inod).*fun(Jnod)...  ;   % t2
                    ).*detJw;


                Kyh(:,Inod,Jnod)=Kyh(:,Inod,Jnod)...
                    +(etaint.*(4*eyy+2*exx).*Deriv(:,2,Inod).*fun(Jnod)...
                    +etaint.*2.*exy.*Deriv(:,1,Inod).*fun(Jnod)...
                    +dtauydh.*funI.*fun(Jnod)...   % +deltaint.*beta2int.*vint.*fun(Inod).*fun(Jnod)...
                    +ca*g*( rhoint - rhow * Dddhint).* dbdy .*funI.*fun(Jnod)...                 % t1:
                    -ca*g*(rhoint.*hint-rhow*dint.*Dddhint).*Deriv(:,2,Inod).*fun(Jnod)...  ;    % t2
                    ).*detJw;

            end


          
            
            Khu(:,Inod,Jnod)=Khu(:,Inod,Jnod)...
                +theta*(rhoint.*dhdx.*fun(Jnod)+drhodx.*hint.*fun(Jnod)+rhoint.*hint.*Deriv(:,1,Jnod))...
                .*SUPG.*detJw*dt; % +dSUPGu.*detJw*dt;
            
            
            Khv(:,Inod,Jnod)=Khv(:,Inod,Jnod)...
                +theta*(rhoint.*dhdy.*fun(Jnod)+drhody.*hint.*fun(Jnod)+rhoint.*hint.*Deriv(:,2,Jnod))...
                .*SUPG.*detJw*dt; % +dSUPGv.*detJw*dt;
    
            Khh(:,Inod,Jnod)=Khh(:,Inod,Jnod)...
                +(rhoint.*fun(Jnod)...
                -dt*theta*rhoint.*dadhint.*fun(Jnod)...
                +dt*theta*rhoint.*fun(Jnod).*h1barr/lambda_h...
                +dt*theta.*(rhoint.*exx.*fun(Jnod)+drhodx.*uint.*fun(Jnod)+rhoint.*uint.*Deriv(:,1,Jnod)+...
                rhoint.*eyy.*fun(Jnod)+drhody.*vint.*fun(Jnod)+rhoint.*vint.*Deriv(:,2,Jnod)))...
                .*SUPG.*detJw;
            
            
        end
    end
    
    % note R=T-F;
    %  dR/dh  dh = -R
    %  dT/dh-dF/dh=-T+F  or dF/dh-dT/dh=T-F
    
    t1=-ca*g*(rhoint.*hint-rhow*dint).*dbdx.*funI+ rhoint.*g.*hint.*sa.*funI;
    % t1=-F.g* Heint.*(rhoint.*hint-F.rhow*dint).*dBdx.*fun(Inod)*ca+ rhoint.*F.g.*hint.*sa.*fun(Inod); % testing
    t2=0.5*ca*g.*(rhoint.*hint.^2-rhow.*dint.^2).*Deriv(:,1,Inod);
    t3=hint.*etaint.*(4*exx+2*eyy).*Deriv(:,1,Inod);
    t4=hint.*etaint.*2.*exy.*Deriv(:,2,Inod);
    
    t5=taux.*funI; % beta2int.*uint.*fun(Inod);
    
    Tx(:,Inod)=Tx(:,Inod)+(t3+t4+t5).*detJw;
    Fx(:,Inod)=Fx(:,Inod)+(t1+t2).*detJw;
    
    t1=-ca*g*(rhoint.*hint-rhow*dint).*dbdy.*funI;
    t2=0.5*ca*g.*(rhoint.*hint.^2-rhow.*dint.^2).*Deriv(:,2,Inod);
    t3=hint.*etaint.*(4*eyy+2*exx).*Deriv(:,2,Inod);
    t4=hint.*etaint.*2.*exy.*Deriv(:,1,Inod);
    t5=tauy.*funI; % beta2int.*vint.*fun(Inod);
    
    Ty(:,Inod)=Ty(:,Inod)+(t3+t4+t5).*detJw;
    Fy(:,Inod)=Fy(:,Inod)+(t1+t2).*detJw;
    
    % Lrhow = 1000/(L*rhow) ; 
    % aw=invLrhow*    dt* (theta*(u0int.*taux+v0int.*tauy)+(1-theta)*(u2int.*taux+v1int.*tauy))

    qterm=  dt*(theta*qx1dx+(1-theta)*qx0dx+theta*qy1dy+(1-theta)*qy0dy).*SUPG;
    dhdt=  rhoint.*(h0int-hint+dt*(1-theta)*h0barr+dt*theta*h1barr).*SUPG;
    accterm=  dt*rhoint.*((1-theta)*a0int+theta*a1int).*SUPG;
    
    th=-dhdt;
    fh=  accterm - qterm;
    
    % R is calculated as R=th-fh  and then I solve K x = -R
    % thus: th has opposite sign but fh not
    % second and third-order Taylor terms
    
    %
    Th(:,Inod)=Th(:,Inod)+th.*detJw;
    Fh(:,Inod)=Fh(:,Inod)+fh.*detJw;
    
    
end
