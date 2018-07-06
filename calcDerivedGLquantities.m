
function [Fx,Fy,Ngl,Tgl,Theta,Ffree,f,qgl,nx,ny,xGLele,yGLele,GLele,hGL,uGL,vGL]=...
    calcDerivedGLquantities(CtrlVar,s,b,h,S,B,u,v,rho,rhow,g,coordinates,connectivity,nip,GF,GLgeo,xGL,yGL,txx,tyy,txy)
    
    % Calculates various quantaties that are of interest for the grounding line
    %
    % All returned variables are element based
    %
    % if CtlVar.calcDerivedGLquantitiesForGLeleOnly==1, then quantities
    % are only calculated for GL elements (as based on GLgeo(:,1)),
    % othewise returned variables defined for the whole domain.
    %
    % The normal to the grounding line is based on the spatial derivatives of the GL.node field (averaged over element))
    %
    % [nx,ny]                            : normals to the GL.node field averaged over elements
    % Fx=(2*txx+tyy).*nx+txy.*ny         :
    % Fy=txy.*nx+(2*tyy+txx).*ny         :
    % Ngl                                : normal component of the [Fx,Fy] vector
    % Tgl                                : tangential component of the [Fx,Fy] vector
    % Theta=Ngl./Ffree                   : Buttressing ratio  (1 indicates no buttressing)
    % Ffree=0.5*g*rho.*(1-rho./rhow).*h  :
    % f=rho.*h./(rhow.*H);               : floating ratio (1 if fully floating)
    % qgl                                : flux normal to GL.node field   (units: kg/m/a )
    % GLele                              : (indices of)) grounding line elements
    
    %
    % the difference between (xGLele,yGLele) and (xGL,yGL) is that the former refers to the element through wich the 
    % grounding line goes, whereas the latter is calculated from the position of the grounding line itself as given by GLgeometry.m
    % (xGL,yGL) is more accurate, but (xGLele,yGLele) reflects better the element-averaged aspect of the quantities calculated.
    %
    % To plot Ngl field for the GL elements only: Ngl(GLgeo(:,1))) or Ngl(GLele)
    %
    
    %%
    %  Force balance at a floating calving front requires
    %  T n = 0.5 \rho (1-\rho/\rho_w)  h  n
    %  where
    %  T=( 2 \txx +  \ tyy          &          \txy         )
    %    (     \txy                 &     2 \tyy + \txx     )
    %  and n is a normal to the grounding line
    %
    %
    % F_{front} =T n = ( 2 \txx + \tyy) nx + \txy ny  )
    %                  ( \txy nx + (2 \tyy + \txx) ny )
    %
    %  At the calving from T n = F_{n} n where Fn is a scalar, i.e. the traction is normal to the front
    %  therefore n' T n = F_{n}= 0.5 \rho (1-\rho/\rho_w)  h
    %  and therefore t' T n = F_{n}= 0  where t is a unit vector in the horizontal plane and tangential to n
    %
    % d=He(h-rhow*(S-B)/rhow)*(S-b)
    %
    %%
    
    [Nele,nod]=size(connectivity);
    
    
    
    % calculationg nx and ny to the grounding line from the gradient of the GL.node field
    [dfdx,dfdy,xint,yint]=calcFEderivatives(GF.node,coordinates,connectivity,nip,CtrlVar);
    %dfdx=repmat(mean(dfdx,2),1,nip); dfdy=repmat(mean(dfdy,2),1,nip);
    temp=sqrt(dfdx.*dfdx+dfdy.*dfdy) ; nx=-dfdx./temp ; ny=-dfdy./temp;
    
    
    nx(~isfinite(nx))=0;  ny(~isfinite(ny))=0;
    
    
    Fx=(2*txx+tyy).*nx+txy.*ny ;  Fy=txy.*nx+(2*tyy+txx).*ny ;     % traction normal to the GF.node field
    
    Ngl=nx.*Fx+ny.*Fy;   Ngl=mean(Ngl,2);                          % normal traction (mean value for each element)
    Tgl=-ny.*Fx+nx.*Fy;   Tgl=mean(Tgl,2);                         % tangential traction (mean value for each element)
    
    xGLele=mean(xint,2);
    yGLele=mean(yint,2);

    
    % normal stresses at a freely floating calving front
    
    
    Ffree=0.5*g*rho.*(1-rho./rhow).*h ;  Ffree=mean(reshape(Ffree(connectivity,1),Nele,nod),2);
    
    
    Theta=Ngl./Ffree;                     % Buttressing ratio, calculated for all elementsfigure ; semilogy(hgl,qgl/1e9,'o') ; xlabel('h at GL (m)') ; ylabel('q at GL (10^9 kg m^{-1} a^{-1})')   ;
    
    % Where fully grounded or normal traction tiny, set to NaN
    %ind=GF.ele>0.999     ; Theta(ind)=NaN; Ffree(ind)=NaN; Ngl(ind)=NaN;
    %ind=abs(Ngl) < eps   ; Theta(ind)=NaN; Ffree(ind)=NaN; Ngl(ind)=NaN;
    
    
    % nodal flux
    qx=u.*h.*rho ; qy=v.*h.*rho ;
    qxEle=mean(reshape(qx(connectivity,1),Nele,nod),2); qyEle=mean(reshape(qy(connectivity,1),Nele,nod),2);
    
    d=S-b;
    f=rho.*h./(rhow*d);    f=mean(reshape(f(connectivity,1),Nele,nod),2);
    
    nxEle=mean(nx,2) ;  nyEle=mean(ny,2) ; 
    qgl=qxEle.*nxEle+qyEle.*nyEle;
    
    nx=mean(nx,2) ;  ny=mean(ny,2) ; 
    GLele=GLgeo(:,1);
            
     if CtrlVar.calcDerivedGLquantitiesForGLeleOnly==1
         Fx=Fx(GLele); Fy=Fy(GLele); Ngl=Ngl(GLele); Tgl=Tgl(GLele); Theta=Theta(GLele); f=f(GLele); 
         nx=nx(GLele) ; ny=ny(GLele); 
     end
     
     qgl=qgl(GLele) ; 
     xGLele=xGLele(GLele) ;  yGLele=yGLele(GLele) ;
     
     x=coordinates(:,1); y=coordinates(:,2); 
     DTxy = DelaunayTri(x,y); 
     hGL=Grid1toGrid2(DTxy,h,xGL,yGL);
     uGL=Grid1toGrid2(DTxy,u,xGL,yGL);
     vGL=Grid1toGrid2(DTxy,v,xGL,yGL);
     
     
    
    return
    
end

