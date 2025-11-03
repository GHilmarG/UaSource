function [kv,rh]=NexthTG3Assemble2DMatrix(dt,h0,u0,v0,du0dt,dv0dt,a0,da0dt,u1,v1,a1,da1dt,du1dt,dv1dt,coordinates,connectivity,Boundary,nip,CtrlVar)
    
    % vectorized, BUT `sparse' is used too often (nod x nod times in total)
    
    %Gamma=CtrlVar.Gamma;
    
    
    Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
    ndim=2; dof=1; neq=dof*Nnodes;
    
    
    h0nod=reshape(h0(connectivity,1),Nele,nod);    % Nele x nod
    u0nod=reshape(u0(connectivity,1),Nele,nod);
    v0nod=reshape(v0(connectivity,1),Nele,nod);
    du0dtnod=reshape(du0dt(connectivity,1),Nele,nod);
    dv0dtnod=reshape(dv0dt(connectivity,1),Nele,nod);
    a0nod=reshape(a0(connectivity,1),Nele,nod);
    da0dtnod=reshape(da0dt(connectivity,1),Nele,nod);
    
    u1nod=reshape(u1(connectivity,1),Nele,nod);
    v1nod=reshape(v1(connectivity,1),Nele,nod);
    du1dtnod=reshape(du1dt(connectivity,1),Nele,nod);
    dv1dtnod=reshape(dv1dt(connectivity,1),Nele,nod);
    a1nod=reshape(a1(connectivity,1),Nele,nod);
    da1dtnod=reshape(da1dt(connectivity,1),Nele,nod);
    
    % can I not just use h0 for h0nod, ie reuse variable?
    
    [points,weights]=sample('triangle',nip,ndim);
    
    kv=sparseUA(neq,neq);
    d1d1=zeros(Nele,nod,nod);
    b1=zeros(Nele,nod);
    
    % vector over all elements for each integartion point
    for Iint=1:nip
        
        
        
        fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        [Deriv,detJ]=derivVector(coordinates,connectivity,nip,Iint);
        
        % Deriv : Nele x dof x nod
        %  detJ : Nele
        
        % values at integration point
        
        h0int=h0nod*fun; u0int=u0nod*fun; v0int=v0nod*fun; du0dtint=du0dtnod*fun;  dv0dtint=dv0dtnod*fun;  a0int=a0nod*fun; da0dtint=da0dtnod*fun;
        u1int=u1nod*fun; v1int=v1nod*fun; du1dtint=du1dtnod*fun;  dv1dtint=dv1dtnod*fun;  a1int=a1nod*fun; da1dtint=da1dtnod*fun;
        
        du1dx=zeros(Nele,1); du0dx=zeros(Nele,1); dh0dx=zeros(Nele,1);
        dv1dy=zeros(Nele,1); dv0dy=zeros(Nele,1); dh0dy=zeros(Nele,1);
        
        % derivatives at integration points
        for Inod=1:nod
            du1dx=du1dx+Deriv(:,1,Inod).*u1nod(:,Inod);
            du0dx=du0dx+Deriv(:,1,Inod).*u0nod(:,Inod);
            
            dv1dy=dv1dy+Deriv(:,2,Inod).*v1nod(:,Inod);
            dv0dy=dv0dy+Deriv(:,2,Inod).*v0nod(:,Inod);
            
            dh0dx=dh0dx+Deriv(:,1,Inod).*h0nod(:,Inod);
            dh0dy=dh0dy+Deriv(:,2,Inod).*h0nod(:,Inod);
            
        end
        
        detJw=detJ*weights(Iint);
        
        
        
        dthalf=dt/2; dt2=dt*dt/12;
        
        for Inod=1:nod
            for Jnod=1:nod
                
                %1st order terms
                h1term=fun(Jnod).*fun(Inod);
                
                h1dxu1=dthalf*du1dx.*fun(Jnod).*fun(Inod);
                u1dxh1=dthalf*u1int.*Deriv(:,1,Jnod).*fun(Inod);
                
                h1dyv1=dthalf*dv1dy.*fun(Jnod).*fun(Inod);
                v1dyh1=dthalf*v1int.*Deriv(:,2,Jnod).*fun(Inod);
                
                %2nd and 3rd order terms
                
                du1dth1=du1dtint.*fun(Jnod).*Deriv(:,1,Inod);
                
                u1u1dh1dx=-u1int.*u1int.*Deriv(:,1,Inod).*Deriv(:,1,Jnod);
                u1du1dxh1=-u1int.*du1dx.*fun(Jnod).*Deriv(:,1,Inod);
                
                u1v1dh1dy=-u1int.*v1int.*Deriv(:,2,Inod).*Deriv(:,1,Jnod);
                u1dv1dyh1=-u1int.*dv1dy.*fun(Jnod).*Deriv(:,1,Inod);
                
                
                
                dv1dth1=dv1dtint.*fun(Jnod).*Deriv(:,2,Inod);
                
                v1u1dh1dx=-v1int.*u1int.*Deriv(:,1,Inod).*Deriv(:,2,Jnod);
                v1du1dxh1=-v1int.*du1dx.*fun(Jnod).*Deriv(:,2,Inod);
                
                v1v1dh1dy=-v1int.*v1int.*Deriv(:,2,Inod).*Deriv(:,2,Jnod);
                v1dv1dyh1=-v1int.*dv1dy.*fun(Jnod).*Deriv(:,2,Inod);
                
                % adding it all up
                d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)...
                    +(h1term+h1dxu1+u1dxh1+h1dyv1+v1dyh1).*detJw...
                    +CtrlVar.TG3*dt2*(du1dth1+u1u1dh1dx+u1du1dxh1+u1v1dh1dy+u1dv1dyh1).*detJw...
                    +CtrlVar.TG3*dt2*(dv1dth1+v1u1dh1dx+v1du1dxh1+v1v1dh1dy+v1dv1dyh1).*detJw;
                
                
            end
            
            % 1st order terms
            h0term=h0int.*fun(Inod);
            aterm=dthalf*(a0int+a1int).*fun(Inod);
            
            
            hdxu0=-dthalf*du0dx.*h0int.*fun(Inod);
            udxh0=-dthalf*dh0dx.*u0int.*fun(Inod);
            
            hdyv0=-dthalf*dv0dy.*h0int.*fun(Inod);
            vdyh0=-dthalf*dh0dy.*v0int.*fun(Inod);
            
            %2nd and 3rd order terms
            
            %x
            du0dth0=du0dtint.*h0int.*Deriv(:,1,Inod);
            
            
            u0u0dh0dx=-u0int.*u0int.*dh0dx.*Deriv(:,1,Jnod);
            u0du0dxh0=-u0int.*du0dx.*h0int.*Deriv(:,1,Inod);
            
            u0v0dh0dy=-u0int.*v0int.*dh0dy.*Deriv(:,1,Jnod);
            u0dv0dyh0=-u0int.*dv0dy.*h0int.*Deriv(:,1,Inod);
            
            u0a0=u0int.*a0int.*Deriv(:,1,Inod);
            u1a1=-u1int.*a1int.*Deriv(:,1,Inod);
            
            
            %y
            dv0dth0=dv0dtint.*h0int.*Deriv(:,2,Inod);
            
            v0u0dh0dx=-v0int.*u0int.*dh0dx.*Deriv(:,2,Jnod);
            v0du0dxh0=-v0int.*du0dx.*h0int.*Deriv(:,2,Inod);
            
            v0v0dh0dy=-v0int.*v0int.*dh0dy.*Deriv(:,2,Jnod);
            v0dv0dyh0=-v0int.*dv0dy.*h0int.*Deriv(:,2,Inod);
            
            v0a0=v0int.*a0int.*Deriv(:,2,Inod);
            v1a1=-v1int.*a1int.*Deriv(:,2,Inod);
            
            % adding it all up
            
            b1(:,Inod)=b1(:,Inod)+(h0term+aterm+hdxu0+udxh0+hdyv0+vdyh0).*detJw...
                +CtrlVar.TG3*dt2*(u0a0+u1a1+v0a0+v1a1).*detJw...
                +CtrlVar.TG3*dt2*(du0dth0+u0u0dh0dx+u0du0dxh0+u0v0dh0dy+u0dv0dyh0).*detJw...
                +CtrlVar.TG3*dt2*(dv0dth0+v0u0dh0dx+v0du0dxh0+v0v0dh0dy+v0dv0dyh0).*detJw...
                +CtrlVar.TG3*dt2*(da0dtint-da1dtint).*detJw;
            
            
        end
    end
    
    % assemble right-hand side
    
    rh=sparseUA(neq,1);
    for Inod=1:nod
        rh=rh+sparseUA(connectivity(:,Inod),ones(Nele,1),b1(:,Inod),neq,1);
    end
    
    
    for Inod=1:nod
        for Jnod=1:nod
            kv=kv+sparseUA(connectivity(:,Inod),connectivity(:,Jnod),d1d1(:,Inod,Jnod),neq,neq);
        end
    end
       
    
    %[kvBoundary,rhBoundary]=BoundaryIntegralPrognosticSemiImplicit(coordinates,connectivity,Boundary,h0,u0,v0,u1,v1,a0,a1,dt,CtrlVar);
    %kv=kv+kvBoundary ; rh=rh+rhBoundary ;
    
    
end