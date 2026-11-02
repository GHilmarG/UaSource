function [wSurf,wSurfInt,wBedInt,wBed]=calcVerticalSurfaceVelocity(rho,rhow,h,S,B,b,ub,vb,as,ab,exx,eyy,xint,yint,coordinates,connectivity,nip,CtrlVar,lGL)
      
    if nargin~=18 ; error(' number of input arguments incorrect ') ; end
    
    % calculates vertical velocity at both integration points and nodal points (vectorized)
    % Note: the projection from interpolation points to nodes is done using simple scattered interpolation
    % this could be improved by finding wNod such that \int wNod Np Nq = \int wInt Nq
    % -> wNod=M\ (\int wInt Nq) 
    % 

    
    [Nele,nod]=size(connectivity); ndim=2;
    [points,weights]=sample('triangle',nip,ndim);
    
    hnod=reshape(h(connectivity,1),Nele,nod);   % Nele x nod
    bnod=reshape(b(connectivity,1),Nele,nod);
    unod=reshape(ub(connectivity,1),Nele,nod);
    vnod=reshape(vb(connectivity,1),Nele,nod);
    Snod=reshape(S(connectivity,1),Nele,nod);
    Bnod=reshape(B(connectivity,1),Nele,nod);
    rhonod=reshape(rho(connectivity,1),Nele,nod);
    asnod=reshape(as(connectivity,1),Nele,nod);
    abnod=reshape(ab(connectivity,1),Nele,nod);
    
    % coox=reshape(coordinates(connectivity,1),Nele,nod);
    % cooy=reshape(coordinates(connectivity,2),Nele,nod);
    
    
    wBedInt=zeros(Nele,nip); wSurfInt=zeros(Nele,nip);
    
    % (2 Nov 2026) transition length of the grounding-line mask: the nodal lGL (optional input, as F.lGL) interpolated to the integration points, otherwise lmin
    kHe=(1/CtrlVar.GroundingLineTransition.lmin) ; lGLnod=[] ;
    if nargin>=19 && numel(lGL)==size(coordinates,1) ; lGLnod=reshape(lGL(connectivity),size(connectivity,1),size(connectivity,2)) ; end
    for Iint=1:nip                           % loop over integration points
        fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        
        hint=hnod*fun;
        uint=unod*fun;
        vint=vnod*fun;
        
        Bint=Bnod*fun;
        Sint=Snod*fun;
        rhoint=rhonod*fun;
        asint=asnod*fun;
        abint=abnod*fun;
        
        hfint=rhow*(Sint-Bint)./rhoint;
        kHi=kHe ; if ~isempty(lGLnod) ; kHi=1./(lGLnod*fun) ; end
        Heint = HeavisideApprox(kHi,hint-hfint,CtrlVar.Hh0);   % (2 Nov 2026) as in the uv assembly
        
        
        [Deriv]=derivVector(coordinates,connectivity,nip,Iint); % Nele x dof x nod
        % The derivative depends on the det of each element
        % is therefore a much bigger array then fun
        
        
        dbdxint=zeros(Nele,1) ; dbdyint=zeros(Nele,1) ;
        for I=1:nod
            dbdxint=dbdxint+Deriv(:,1,I).*bnod(:,I);
            dbdyint=dbdyint+Deriv(:,2,I).*bnod(:,I);
        end
        
        wBedInt(:,Iint)=uint.*dbdxint+vint.*dbdyint;
        %wDiffInt=-hint.*(exx(:,Iint)+eyy(:,Iint)); % difference in surface and basal vertical velocity
        
        % Grounded part:  wSurf=ab-h (exx+eyy) +u \p_x b + v \p_y b  + \p_t B
        wsGround=abint+uint.*dbdxint+vint.*dbdyint-hint.*(exx(:,Iint)+eyy(:,Iint));
        
        % Floating part: wSurf=-(1-rho/rhow) ( h (exx+eyy) - a_b ) -a_s rho/rhow
        wsFloat=-(1-rhoint/rhow).*(hint.*(exx(:,Iint)+eyy(:,Iint))-abint)-asint.*rhoint/rhow;
        
        wSurfInt(:,Iint)=wsGround.*Heint+wsFloat.*(1-Heint) ;
        
    end
    
    % now simply interpolate onto nodal points
    Xint=xint(:) ; Yint=yint(:);
    [~, Iint, ~] = unique([Xint Yint],'first','rows');
    Iint = sort(Iint); Xint = Xint(Iint); Yint = Yint(Iint);
    DTint = DelaunayTri(Xint,Yint);
    %ic=incenters(DTint); [cn,on] = inpoly(ic,MeshBoundaryCoordinates); DTintTriInside=DTint.Triangulation(cn,:);
    wtemp=wSurfInt(:) ; wtemp=wtemp(Iint);
    wSurf=Grid1toGrid2(DTint,wtemp,coordinates(:,1),coordinates(:,2));
    
    if nargout>3
        wtemp=wBedInt(:) ; wtemp=wtemp(Iint);
        wBed=Grid1toGrid2(DTint,wtemp,coordinates(:,1),coordinates(:,2));
    end
    
    
end
