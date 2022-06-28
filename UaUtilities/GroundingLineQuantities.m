function GLQ=GroundingLineQuantities(F,xGL,yGL,nxGL,nyGL)


% calculates various quantities of interest along GL
% This is done by:
% 1) defining a grounding line curve using splines
% 2) interpolating nodal and integration-point fields onto this curve
% 3) calculating other quantities of interest

% for each grounding line do:


%    GLQ.ThetaN=GLQ.N./GLQ.Ffree;                     % (normal) Buttressing fraction
%    GLQ.ThetaT=GLQ.T./GLQ.Ffree;                     % (tangential) Buttressing fraction



    % interpolate onto GL curve
    
    Interpolant=scatteredInterpolant;

    Interpolant.Points=[F.x F.y] ;

    Interpolant.Values=F.ub ; ugl=Interpolant(xGL,yGL) ; 
    Interpolant.Values=F.vb ; vgl=Interpolant(xGL,yGL) ; 
    Interpolant.Values=F.h ; hgl=Interpolant(xGL,yGL) ; 
    Interpolant.Values=F.rho ; rhogl=Interpolant(xGL,yGL) ; 

    Interpolant.Values=F.txx ; txxgl=Interpolant(xGL,yGL) ; 
    Interpolant.Values=F.txy ; txygl=Interpolant(xGL,yGL) ; 
    Interpolant.Values=F.tyy ; tyygl=Interpolant(xGL,yGL) ; 

    % d=S-b;
    % [hgl,ugl,vgl,rhogl,dgl]=InterpolateNodalVariables(DTxy,xGL,yGL,h,u,v,rho,d);
    % [txxgl,tyygl,txygl]=InterpolateIntVariables(DTint,Iint,xGL,yGL,txx,tyy,txy);
    
    
    GLQ.Fx=(2*txxgl+tyygl).*nxGL+txygl.*nyGL ;
    GLQ.Fy=txygl.*nxGL+(2*tyygl+txxgl).*nyGL ;
    
    GLQ.N=nxGL.*GLQ.Fx+nyGL.*GLQ.Fy;           % normal traction
    GLQ.T=-nyGL.*GLQ.Fx+nxGL.*GLQ.Fy;          % tangential traction
    
    % normal stresses at a freely floating calving front
    
    
     GLQ.Ffree=0.5*F.g*rhogl.*(1-rhogl./F.rhow).*hgl ;  % Note that here 0.5 is correct
    % in 1D we get 
    %
    % (2*txxgl+tyygl).*nxGL+txygl.*nyGL = (2*txxgl+0).*1+0.*0 ;
    %                                   =  2*txxgl
    % 
    % or   txxgl = 0.5 g rho (1-rho/rhow) h /2
    %            = 0.25 g rho (1-rho/rhow) h 
    
    
    
    GLQ.ThetaN=GLQ.N./GLQ.Ffree;                     % (normal) Buttressing fraction
    GLQ.ThetaT=GLQ.T./GLQ.Ffree;                     % (tangential) Buttressing fraction
    
    % Theta_n is normal component of the horizontal traction, normalized by the 1D value
    % 
    % Theta_n=1 implies no butressing, 
    %
    % Theta>1 implies `positive buttressing' ie the normal (tensional) stress
    % at the grounding line is less than the (tensional) stress in the case of an unconfined ice shelf,
    % in other words the ice flow is restricted by the ice shelf

    % Theta<1 implies normal (tensional) stress at the grounding line that is larger than for the unconfined ice shelf case, ie the
    % ice is pulled out
    
    
    GLQ.kappaN=(GLQ.Ffree-GLQ.N)./GLQ.Ffree;                     % (normal) Buttressing number
    GLQ.kappaT=GLQ.T./GLQ.Ffree;                                 % (tangential) Buttressing number
    
    GLQ.qn=rhogl.*hgl.*(ugl.*nxGL+vgl.*nyGL);        % normal flux
    GLQ.qt=rhogl.*hgl.*(-ugl.*nyGL+vgl.*nxGL);       % tangential flux
    
    
    
    % f=(rhogl.*hgl-rhow*dgl)./(rhogl.*hgl); % floating ratio, 0 fully floating, 1 grounded
    
    
    
end

