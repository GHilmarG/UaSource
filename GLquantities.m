function [qngl,qtgl,ThetaN,ThetaT,kappaN,kappaT,Ngl,Tgl,f,hgl,dgl,ugl,vgl,rhogl,txxgl,tyygl,txygl]=...
        GLquantities(CtrlVar,xGL,yGL,nxGL,nyGL,coordinates,DTxy,DTint,Iint,s,b,S,B,h,u,v,g,rho,rhow,txx,tyy,txy)
    
    
    % calculates various quantities of interest along GL
    % This is done by:
    % 1) defining a grounding line curve using splines
    % 2) interpolating nodal and integration-point fields onto this curve
    % 3) calculating other quantities of interest
    
    % for each grounding line do:
    
    % Note: Old and outdated, consider using GroundingLineQuantities.m instead
    
    
    
    % interpolate onto GL curve
    
    d=S-b;
    [hgl,ugl,vgl,rhogl,dgl]=InterpolateNodalVariables(DTxy,xGL,yGL,h,u,v,rho,d);
    [txxgl,tyygl,txygl]=InterpolateIntVariables(DTint,Iint,xGL,yGL,txx,tyy,txy);
    
    
    Fx=(2*txxgl+tyygl).*nxGL+txygl.*nyGL ;
    Fy=txygl.*nxGL+(2*tyygl+txxgl).*nyGL ;
    
    Ngl=nxGL.*Fx+nyGL.*Fy;           % normal traction
    Tgl=-nyGL.*Fx+nxGL.*Fy;          % tangential traction
    
    % normal stresses at a freely floating calving front
    
    
    Ffree=0.5*g*rhogl.*(1-rhogl./rhow).*hgl ;
    
    
    ThetaN=Ngl./Ffree;                     % (normal) Buttressing fraction
    ThetaT=Tgl./Ffree;                     % (tangential) Buttressing fraction
    
    % Theta_n is difference between the actual normal stress at the grounding line
    % and it's 1d value for an unconfined ice shelf, normalized by the 1d value
    % Theta_n=0 implies no butressing, 
    % Theta>0 implies `positive buttressing' ie the normal (tensional) stress
    % at the grounding line is less than the (tensional) stress in the case of an unconfined ice shelf,
    % in other words the ice flow is restricted by the ice shelf
    % Theta<0 implies normal (tensional) stress at the grounding line that is larger than for the unconfined ice shelf case, ie the
    % ice is pulled out
    
    
    kappaN=(Ffree-Ngl)./Ffree;                     % (normal) Buttressing number
    kappaT=Tgl./Ffree;                             % (tangential) Buttressing number
    
    qngl=rhogl.*hgl.*(ugl.*nxGL+vgl.*nyGL);        % normal flux
    qtgl=rhogl.*hgl.*(-ugl.*nyGL+vgl.*nxGL);       % tangential flux
    
    
    
    f=(rhogl.*hgl-rhow*dgl)./(rhogl.*hgl); % floating ratio, 0 fully floating, 1 grounded
    
    
    
end

