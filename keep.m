function [UserVar,rh,kv,Tv,Lv,Pv,Qx,Qy]=LevelSetEquationAssemblyNR2consistent(UserVar,CtrlVar,MUA,f0,c0,u0,v0,f1,c1,u1,v1)
    
  
    narginchk(11,11)
    
    ndim=2; dof=1; neq=dof*MUA.Nnodes;
    
    theta=CtrlVar.LevelSetTheta;
    dt=CtrlVar.dt;
    CtrlVar.Tracer.SUPG.tau=CtrlVar.LevelSetSUPGtau;

    isL=CtrlVar.LSF.L ; isP=CtrlVar.LSF.P ; isT=CtrlVar.LSF.T ;

    f0nod=reshape(f0(MUA.connectivity,1),MUA.Nele,MUA.nod);
    f1nod=reshape(f1(MUA.connectivity,1),MUA.Nele,MUA.nod);
    
    u0nod=reshape(u0(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
    v0nod=reshape(v0(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
    
    u1nod=reshape(u1(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
    v1nod=reshape(v1(MUA.connectivity,1),MUA.Nele,MUA.nod);   % MUA.Nele x nod
    
    
    c0nod=reshape(c0(MUA.connectivity,1),MUA.Nele,MUA.nod);
    c1nod=reshape(c1(MUA.connectivity,1),MUA.Nele,MUA.nod);
    
    d1d1=zeros(MUA.Nele,MUA.nod,MUA.nod);
    %b1=zeros(MUA.Nele,MUA.nod);
    T=zeros(MUA.Nele,MUA.nod);
    P=zeros(MUA.Nele,MUA.nod);
    L=zeros(MUA.Nele,MUA.nod);
    qx=zeros(MUA.Nele,MUA.nod);
    qy=zeros(MUA.Nele,MUA.nod);
    
    
    if CtrlVar.LevelSetSolutionMethod=="Newton Raphson"
        NR=1;
    else
        NR=0;
    end
    
    
    
    
    
    
    
    
    % vector over all elements for each  integration point
    for Iint=1:MUA.nip  %Integration points
        
        
        
        fun=shape_fun(Iint,ndim,MUA.nod,MUA.points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        
        
        Deriv=MUA.Deriv(:,:,:,Iint);
        detJ=MUA.DetJ(:,Iint);
        
  
        f0int=f0nod*fun;
        f1int=f1nod*fun;
        
        u0int=u0nod*fun; v0int=v0nod*fun;
        u1int=u1nod*fun; v1int=v1nod*fun;
        c0int=c0nod*fun;
        c1int=c1nod*fun;
        
   
        
        % derivatives at one integration point for all elements
        df0dx=zeros(MUA.Nele,1); df0dy=zeros(MUA.Nele,1);
        df1dx=zeros(MUA.Nele,1); df1dy=zeros(MUA.Nele,1);
        for Inod=1:MUA.nod
            
            
            df0dx=df0dx+Deriv(:,1,Inod).*f0nod(:,Inod);
            df0dy=df0dy+Deriv(:,2,Inod).*f0nod(:,Inod);
            
            df1dx=df1dx+Deriv(:,1,Inod).*f1nod(:,Inod);
            df1dy=df1dy+Deriv(:,2,Inod).*f1nod(:,Inod);
        end
        
        
        
        
        
        % Norm of gradient (NG)
        NG0=sqrt(df0dx.*df0dx+df0dy.*df0dy); % at each integration point for all elements
        NG1=sqrt(df1dx.*df1dx+df1dy.*df1dy); % at each integration point for all elements
        n1x=-df1dx./NG1;  n1y=-df1dy./NG1;
        n0x=-df0dx./NG0;  n0y=-df0dy./NG0;
        
      % if gradient is very small, set normal to zero, and with it the cx and cy components
        I0=NG0< eps^2 ;
        I1=NG1< eps^2 ;
        n1x(I1)=0 ; n1y(I1)=0;
        n0x(I0)=0 ; n0y(I0)=0;

        
        cx1int=-c1int.*n1x ; cy1int=-c1int.*n1y;
        cx0int=-c0int.*n0x ; cy0int=-c0int.*n0y;
        
        %% limit cx-u and cy-v where it is suffiently far away from the zero level
        
   
        
        %%
        
        tauSUPGint=CalcSUPGtau(CtrlVar,MUA,u0int-cx0int,v0int-cy0int,dt); 
        %tauSUPGint=CalcSUPGtau(CtrlVar,MUA,u0int,v0int,dt); 
        
        % I need to think about a good def for mu
        %
        % Idea :  sqrt( (u0int-cx0int).^2+(v0int-cy0int).^2)) .*sqrt(2*MUA.EleAreas) ;
        %
        
        
        switch CtrlVar.LevelSetFABmu.Scale
            
            case "constant"
                Scale=1 ;
            case "ucl"
                Scale =  sqrt( (u0int-cx0int).^2+(v0int-cy0int).^2) .*sqrt(2*MUA.EleAreas) ;
            otherwise
                
                error(adsf')
        end
        
        mu=Scale*CtrlVar.LevelSetFABmu.Value;  % This had the dimention l^2/t
        
        
        [kappaint0]=LevelSetEquationFAB(CtrlVar,NG0,mu);
        [kappaint1,dkappa]=LevelSetEquationFAB(CtrlVar,NG1,mu);
   
        
        detJw=detJ*MUA.weights(Iint);
        
        if any(~isfinite(n0x)) ||any(~isfinite(n1x)) || any(~isfinite(kappaint1)) || any(~isfinite(dkappa))
            save TestSaveLSFnotFinite
            error("LevelSetEquationAssemblyNR2:notfinite","n0x, n1x kappa not finite")
        end

        
        for Inod=1:MUA.nod
            
          
            %SUPG=fun(Inod)+CtrlVar.Tracer.SUPG.Use*tauSUPGint.*(u0int.*Deriv(:,1,Inod)+v0int.*Deriv(:,2,Inod));
            SUPG=fun(Inod)+CtrlVar.Tracer.SUPG.Use*tauSUPGint.*((u0int-cx0int).*Deriv(:,1,Inod)+(v0int-cy0int).*Deriv(:,2,Inod));
            SUPGdetJw=SUPG.*detJw;
            
            if nargout>2
                for Jnod=1:MUA.nod
                    
                    
                    Tlhs=fun(Jnod).*SUPGdetJw ; 
                        
                    Llhs=...
                        +NR*dt*theta*(...
                        (u1int-cx1int).*Deriv(:,1,Jnod) + (v1int-cy1int).*Deriv(:,2,Jnod))... 
                        .*SUPGdetJw;

                     
                    % Pertubation term
                    Plhs=dt*theta*...
                        +(kappaint1.*(Deriv(:,1,Jnod).*Deriv(:,1,Inod)+Deriv(:,2,Jnod).*Deriv(:,2,Inod)) ...
                        - NR*dkappa.*(n1x.*Deriv(:,1,Jnod)+n1y.*Deriv(:,2,Jnod)).*(df1dx.*Deriv(:,1,Inod)+df1dy.*Deriv(:,2,Inod))) ...
                        .*detJw;
               
           
                    d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+isL*Llhs+isP*Plhs+isT*Tlhs;
                    
                end
            end

            % Note, I solve: 
            %         
            %   LSH  \phi  = - RHS
            %
            
            Trhs=(f1int-f0int).*SUPGdetJw ; 

            Lrhs= ...
                +    dt*theta* ((u1int-cx1int).*df1dx +(v1int-cy1int).*df1dy).*SUPGdetJw...
                + dt*(1-theta)*((u0int-cx0int).*df0dx +(v0int-cy0int).*df0dy).*SUPGdetJw;

 
            
            % Pertubation term
            Prhs=...
                       dt*theta*kappaint1.*(df1dx.*Deriv(:,1,Inod)+df1dy.*Deriv(:,2,Inod)).*detJw ...
                 + dt*(1-theta)*kappaint0.*(df0dx.*Deriv(:,1,Inod)+df0dy.*Deriv(:,2,Inod)).*detJw;
            
            % qx= kappaint0.*df0dx ; 
            % qu= kappaint0.*df0dy) ; 
            qx(:,Inod)=qx(:,Inod)+kappaint0.*df0dx ; 
            qy(:,Inod)=qy(:,Inod)+kappaint0.*df0dy ; 
             
            P(:,Inod)=P(:,Inod)+Prhs; 
            L(:,Inod)=L(:,Inod)+Lrhs; 
            T(:,Inod)=T(:,Inod)+Trhs; 
            
        end
    end
    
    Pv=sparseUA(neq,1);
    Lv=sparseUA(neq,1);
    Tv=sparseUA(neq,1);
    Qx=sparseUA(neq,1);
    Qy=sparseUA(neq,1);
    for Inod=1:MUA.nod
        Pv=Pv+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),P(:,Inod),neq,1);
        Lv=Lv+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),L(:,Inod),neq,1);
        Tv=Tv+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),T(:,Inod),neq,1);
        Qx=Qx+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),qx(:,Inod),neq,1);
        Qy=Qy+sparseUA(MUA.connectivity(:,Inod),ones(MUA.Nele,1),qy(:,Inod),neq,1);
    end
    
    rh=isL*Lv+isP*Pv+isT*Tv; 
    
    if nargout>2
        Iind=zeros(MUA.nod*MUA.nod*MUA.Nele,1); Jind=zeros(MUA.nod*MUA.nod*MUA.Nele,1);Xval=zeros(MUA.nod*MUA.nod*MUA.Nele,1);
        istak=0;
        
        for Inod=1:MUA.nod
            for Jnod=1:MUA.nod
                Iind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Inod);
                Jind(istak+1:istak+MUA.Nele)=MUA.connectivity(:,Jnod);
                Xval(istak+1:istak+MUA.Nele)=d1d1(:,Inod,Jnod);
                istak=istak+MUA.Nele;
            end
        end
        
        kv=sparseUA(Iind,Jind,Xval,neq,neq);
    end
    
end