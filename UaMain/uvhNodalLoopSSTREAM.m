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

%%
%
% Highly vectorized over both the integration points and the first nodal loop.
%
%
%%

Nele=size(Deriv,1);

Dx=reshape(Deriv(:,1,:),Nele,nod);
Dy=reshape(Deriv(:,2,:),Nele,nod);

funR=fun.';                      % 1 x nod

qx1dx=rhoint.*exx.*hint+rhoint.*uint.*dhdx+drhodx.*uint.*hint;
qy1dy=rhoint.*eyy.*hint+rhoint.*vint.*dhdy+drhody.*vint.*hint;

qx0dx=rhoint.*exx0.*h0int+rhoint.*u0int.*dh0dx+drhodx.*u0int.*h0int;
qy0dy=rhoint.*eyy0.*h0int+rhoint.*v0int.*dh0dy+drhody.*v0int.*h0int;

% i-dependent, Nele x nod.  Always needed: used by the residual as well.
SUPG=funR+0.5.*tau0.*(u0int.*Dx+v0int.*Dy);

if ~Ronly

    he=hint.*etaint;
    heW=he.*detJw;

    a1=hint.*(4.*exx+2.*eyy);
    a2=2*hint.*exy;
    a3=hint.*(4.*eyy+2.*exx);

    % j-dependent, Nele x nod
    Deu=Eint.*((2*exx+eyy).*Dx+exy.*Dy);
    Dev=Eint.*((2*eyy+exx).*Dy+exy.*Dx);

    % i-dependent, Nele x nod, quadrature weight folded in
    GuW=(a1.*Dx+a2.*Dy).*detJw;
    GvW=(a3.*Dy+a2.*Dx).*detJw;

    tuuW=dtauxdu.*detJw;
    tvvW=dtauydv.*detJw;
    tuvW=dtauxdv.*detJw;
    tvuW=dtauydu.*detJw;

    % --- the two momentum-vs-h blocks, rank one: Kxh = Px_i fun_j --------

    t2coeff=ca*g*(rhoint.*hint-rhow*dint.*Dddhint);      % multiplies dphid{x,y}_i

    % if CtrlVar.Development.Pre2025uvhAssembly  % it will be possible to get rid of this...
    %     fx=dtauxdh ...
    %         + ca*g*rhoint.*Heint.*dBdx ...
    %         + ca*g*deltaint.*(rhoint.*hint-rhow*Hposint).*dBdx ...
    %         - sa*g*rhoint;
    %     fy=dtauydh ...
    %         + ca*g*rhoint.*Heint.*dBdy ...
    %         + ca*g*deltaint.*(rhoint.*hint-rhow*Hposint).*dBdy;
    % else
    fx=dtauxdh + ca*g*(rhoint-rhow*Dddhint).*dbdx - sa*g*rhoint;
    fy=dtauydh + ca*g*(rhoint-rhow*Dddhint).*dbdy;
    %end

    PxW=((etaint.*(4*exx+2*eyy)-t2coeff).*Dx + etaint.*2.*exy.*Dy + fx.*funR).*detJw;
    PyW=((etaint.*(4*eyy+2*exx)-t2coeff).*Dy + etaint.*2.*exy.*Dx + fy.*funR).*detJw;

    % --- the three thickness rows, rank one: Kh* = SUPG_i Q*_j ----------

    QuW=((rhoint.*dhdx+drhodx.*hint).*funR + rhoint.*hint.*Dx)*theta.*detJw*dt;
    QvW=((rhoint.*dhdy+drhody.*hint).*funR + rhoint.*hint.*Dy)*theta.*detJw*dt;

    QhW=( (rhoint - dt*theta*rhoint.*dadhint + dt*theta*rhoint.*h1barr/lambda_h ...
           + dt*theta.*(rhoint.*exx+drhodx.*uint+rhoint.*eyy+drhody.*vint)).*funR ...
          + dt*theta.*rhoint.*uint.*Dx + dt*theta.*rhoint.*vint.*Dy ).*detJw;

    for Inod=1:nod

        dxi=Dx(:,Inod);  dyi=Dy(:,Inod);
        gu=GuW(:,Inod);  gv=GvW(:,Inod);
        supi=SUPG(:,Inod);
        PP=fun(Inod)*funR;

        Kxu(:,:,Inod)=Kxu(:,:,Inod)+ 4*heW.*dxi.*Dx + heW.*dyi.*Dy + gu.*Deu + tuuW.*PP;
        Kyv(:,:,Inod)=Kyv(:,:,Inod)+ 4*heW.*dyi.*Dy + heW.*dxi.*Dx + gv.*Dev + tvvW.*PP;
        Kxv(:,:,Inod)=Kxv(:,:,Inod)+ heW.*(2*dxi.*Dy+dyi.*Dx)      + gu.*Dev + tuvW.*PP;
        Kyu(:,:,Inod)=Kyu(:,:,Inod)+ heW.*(2*dyi.*Dx+dxi.*Dy)      + gv.*Deu + tvuW.*PP;

        Kxh(:,:,Inod)=Kxh(:,:,Inod)+ PxW(:,Inod).*funR;
        Kyh(:,:,Inod)=Kyh(:,:,Inod)+ PyW(:,Inod).*funR;

        Khu(:,:,Inod)=Khu(:,:,Inod)+ supi.*QuW;
        Khv(:,:,Inod)=Khv(:,:,Inod)+ supi.*QvW;
        Khh(:,:,Inod)=Khh(:,:,Inod)+ supi.*QhW;

    end

end

%% residual and external forces, vectorised over Inod

c1x=(-ca*g*(rhoint.*hint-rhow*dint).*dbdx + rhoint.*g.*hint.*sa).*detJw;
c1y=(-ca*g*(rhoint.*hint-rhow*dint).*dbdy).*detJw;
c2 =(0.5*ca*g.*(rhoint.*hint.^2-rhow.*dint.^2)).*detJw;
c3x=(hint.*etaint.*(4*exx+2*eyy)).*detJw;
c3y=(hint.*etaint.*(4*eyy+2*exx)).*detJw;
c4 =(hint.*etaint.*2.*exy).*detJw;
txW=taux.*detJw;
tyW=tauy.*detJw;

Tx=Tx + c3x.*Dx + c4.*Dy + txW.*funR;
Fx=Fx + c1x.*funR + c2.*Dx;

Ty=Ty + c3y.*Dy + c4.*Dx + tyW.*funR;
Fy=Fy + c1y.*funR + c2.*Dy;

qterm  = dt*(theta*qx1dx+(1-theta)*qx0dx+theta*qy1dy+(1-theta)*qy0dy);
dhdtq  = rhoint.*(h0int-hint+dt*(1-theta)*h0barr+dt*theta*h1barr);
accterm= dt*rhoint.*((1-theta)*a0int+theta*a1int);

Th=Th + (-dhdtq.*detJw).*SUPG;
Fh=Fh + ((accterm-qterm).*detJw).*SUPG;

end

