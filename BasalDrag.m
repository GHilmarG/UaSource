function [taubx,tauby,dtaubxdu,dtaubxdv,dtaubydu,dtaubydv,dtaubxdh,dtaubydh,taubxo,taubyo,taubxa,taubya] = ...
        BasalDrag(CtrlVar,MUA,He,delta,h,B,H,rho,rhow,ub,vb,C,m,uo,vo,Co,mo,ua,va,Ca,ma,q,g,muk)
    
    narginchk(24,24)
    
    %%
    % Returns basal drag and the directional derivatives of basal drag with respect to u,
    % v and h.
    %
    % taubx     :  x component of basal drag
    % tauby     :  y component of basal drag
    % dtaubxdu  :  derivative of x component of basal drag with respect to u
    % dtaubxdv  :  derivative of x component of basal drag with respect to v
    % dtaubxdh  :  derivative of x component of basal drag with respect to h
    % dtaubydu  :  derivative of y component of basal drag with respect to u
    % dtaubydv  :  derivative of y component of basal drag with respect to v
    % dtaubydh  :  derivative of y component of basal drag with respect to h
    %
    %   [ Delta tbx ] =  M    [Delta u]
    %   [ Delta tby ]         [Delta v]
    %                         [Delta h]
    %
    % where
    %
    %  M = [ dtaubxdu   dtaubxdv   dtaubxdh ]
    %      [ dtaubydu   dtaubydv   dtaubydh ]
    %
    %  Regularisation parameters:
    %
    %  CtrlVar.Czero
    %  CtrlVar.SpeedZero
    %  CtrlVar.HeZero
    %
    %
    %
    %      He = HeavisideApprox(kH,h-hf,CtrlVar.Hh0);
    %   delta = DiracDelta(kH,h-hf,CtrlVar.Hh0);
    %      hf = rhow*H./rho;
    %       H = S-B
    %
    %%
    
    Nouts=nargout;
    
    
    %% Basal drag term : ice
    % this drag term is zero if the velocities are zero.
    
    
    U=(sqrt(ub.*ub+vb.*vb+CtrlVar.SpeedZero^2)).^(1./m-1) ;
    beta2i=(C+CtrlVar.Czero).^(-1./m).*U ; %   (sqrt(ub.*ub+vb.*vb+CtrlVar.SpeedZero^2)).^(1./m-1) ;
    
    %
    %
    %
    
    
    % Dbeta2i is zero for m=1.
    Dbeta2i=(1./m-1).*(C+CtrlVar.Czero).^(-1./m).*(ub.^2+vb.^2+CtrlVar.SpeedZero^2).^((1-3*m)./(2*m));
    He=He+CtrlVar.HeZero ; % Regularisation
    
    
    if ~isfield(CtrlVar,"SlidingLaw")
        CtrlVar.SlidingLaw="tauPower";
    end
    
    
    if CtrlVar.Inverse.dFuvdClambda
        
        switch CtrlVar.SlidingLaw
            case {"Weertman","tauPower"}
                
                % tau = He * (C+CtrlVar.Czero).^(-1./m)   * U
                %  just take the derivative with respect to C
                %
                dFuvdC =  He.*    (1./m).*(C+CtrlVar.Czero).^(-1./m-1)   .*U;
                
            case {"Budd","tauPowerNperfectPower"}
                
                
                %             % tau = Nqm * (C+CtrlVar.Czero).^(-1./m)   * U
                %
                %             hf=rhow.*H./rho;
                %             Dh=h-hf; Dh(Dh<eps)=0;
                %             N=He.*rho.*g.*Dh ;
                %             qm=q./m;
                %             Nqm=N.^(qm) ;
                %
                %             dFuvdC= Nqm.*(1./m).*(C+CtrlVar.Czero).^(-1./m-1)  .*U;
                
                
                %  just take the derivative with respect to C
                %   tau = He.*Nqm.* beta2i
                %       = He.*Nqm.* (C+CtrlVar.Czero).^(-1./m).*U ; %   with U=(sqrt(ub.*ub+vb.*vb+CtrlVar.SpeedZero^2)).^(1./m-1) ;
                %
                
                
                
                hf=rhow.*H./rho;
                hf(hf<eps)=0;
                Dh=h-hf; Dh(Dh<eps)=0;
                N=rho.*g.*Dh;
                qm=q./m;
                Nqm=N.^(qm) ;
                
                
                dFuvdC= He.*Nqm  .*(1./m).*(C+CtrlVar.Czero).^(-1./m-1)  .*U;
                
                
            case {"Tsai"}
                
                
                fprintf("Inversion using Tsai sliding law not implemented. \n")
                error("BasalDrag:InvalidCase","Inversion using Tsai sliding law not implemented.")
                
            otherwise
                
                error("BasalDrag:InvalidCase","Unknown case.")
                
        end
        taubx=dFuvdC;
        return
    end
    
    
    switch CtrlVar.SlidingLaw
        
        case {"Weertman","tauPower"}
            
            [taubxi,taubyi,dtaubxdui,dtaubxdvi,dtaubydui,dtaubydvi,dtaubxdhi,dtaubydhi] = Weertman() ;
            
            
        case {"Budd","tauPowerNperfectPower"}
            
            [N,dNdh]=NPerfect();
            [taubxi,taubyi,dtaubxdui,dtaubxdvi,dtaubydui,dtaubydvi,dtaubxdhi,dtaubydhi] = Budd() ;

            
        case "Coulomb"
            
            [N,dNdh]=NPerfect();
            [taubxi,taubyi,dtaubxdui,dtaubxdvi,dtaubydui,dtaubydvi,dtaubxdhi,dtaubydhi] = Coulomb() ;
            
        case {"Tsai"}
            
            
            N=NPerfect();
            [taubxi,taubyi,dtaubxdui,dtaubxdvi,dtaubydui,dtaubydvi,dtaubxdhi,dtaubydhi] =Weertman();
            [taubxiC,taubyiC,dtaubxduiC,dtaubxdviC,dtaubyduiC,dtaubydviC,dtaubxdhiC,dtaubydhiC] = Coulomb();
            
            TauWeertman2=taubxi.*taubxi+taubyi.*taubyi ;
            TauCoulomb2=taubxiC.*taubxiC+taubyiC.*taubyiC ;
            
            % decide where to use Weertman and where Coulomb
            isCoulomb=TauCoulomb2 <  TauWeertman2 ;
            
            % and now replace where needed
            taubxi(isCoulomb)=taubxiC(isCoulomb) ;
            taubyi(isCoulomb)=taubyiC(isCoulomb) ;
            dtaubxdui(isCoulomb)=dtaubxduiC(isCoulomb);
            dtaubydvi(isCoulomb)=dtaubydviC(isCoulomb);
            dtaubxdvi(isCoulomb)=dtaubxdviC(isCoulomb);
            dtaubydui=dtaubxdvi;  % just symmetry, always true, both Weertman and Coulom
            dtaubxdhi(isCoulomb)=dtaubxdhiC(isCoulomb);
            dtaubydhi(isCoulomb)=dtaubydhiC(isCoulomb);
            
            
        otherwise
            
            error("BasalDrag:CaseNotFound","what sliding law?")
    end
    
    
    %% Sea ice drag term : ocean
    
    if CtrlVar.IncludeMelangeModelPhysics
        
        
        U=ub-uo;
        V=vb-vo;
        
        
        beta2o=(Co+CtrlVar.Czero).^(-1./mo).*(sqrt(U.*U+V.*V+CtrlVar.SpeedZero^2)).^(1./mo-1) ;
        
        taubxo=(1-He).*beta2o.*U;
        taubyo=(1-He).*beta2o.*V;
        
        Dbeta2o=(1./mo-1).*(Co+CtrlVar.Czero).^(-1./mo).*(U.^2+V.^2+CtrlVar.SpeedZero^2).^((1-3*mo)./(2*mo));
        
        dtaubxduo=(1-He).*(beta2o+Dbeta2o.*U.*U);
        dtaubxdvo=(1-He).*(Dbeta2o.*U.*V);
        
        dtaubyduo=dtaubxdvo;
        dtaubydvo=(1-He).*(beta2o+Dbeta2o.*V.*V);
        
        dtaubxdho=-delta.*beta2o.*U ;
        dtaubydho=-delta.*beta2o.*V ;
        
        
        
        
        U=ub-ua;
        V=vb-va;
        
        
        beta2a=(Ca+CtrlVar.Czero).^(-1./ma).*(sqrt(U.*U+V.*V+CtrlVar.SpeedZero^2)).^(1./ma-1) ;
        
        taubxa=(1-He).*beta2a.*U;
        taubya=(1-He).*beta2a.*V;
        
        Dbeta2a=(1./ma-1).*(Ca+CtrlVar.Czero).^(-1./ma).*(U.^2+V.^2+CtrlVar.SpeedZero^2).^((1-3*ma)./(2*ma));
        
        dtaubxdua=(1-He).*(beta2a+Dbeta2a.*U.*U);
        dtaubxdva=(1-He).*(Dbeta2a.*U.*V);
        
        dtaubydua=dtaubxdva;
        dtaubydva=(1-He).*(beta2a+Dbeta2a.*V.*V);
        
        dtaubxdha=-delta.*beta2a.*U ;
        dtaubydha=-delta.*beta2a.*V ;
        
        
    else
        
        taubxo=0;
        taubyo=0;
        
        dtaubxduo=0;
        dtaubxdvo=0;
        
        dtaubyduo=0;
        dtaubydvo=0;
        
        dtaubxdho=0;
        dtaubydho=0;
        
        taubxa=0;
        taubya=0;
        
        dtaubxdua=0;
        dtaubxdva=0;
        
        dtaubydua=0;
        dtaubydva=0;
        
        dtaubxdha=0;
        dtaubydha=0;
    end
    
    
    %% Add up
    
    taubx=taubxi+taubxo+taubxa ;
    tauby=taubyi+taubyo+taubya ;
    
    if nargout>2
        dtaubxdu=dtaubxdui+dtaubxduo+dtaubxdua;
        dtaubxdv=dtaubxdvi+dtaubxdvo+dtaubxdva;
        
        dtaubydu=dtaubydui+dtaubyduo+dtaubydua;
        dtaubydv=dtaubydvi+dtaubydvo+dtaubydva;
        
        dtaubxdh=dtaubxdhi+dtaubxdho+dtaubxdha;
        dtaubydh=dtaubydhi+dtaubydho+dtaubydha;
        
    end
    
    %% Nested functions
    function [N,dNdh]=NPerfect()
        
        hf=rhow.*H./rho;
        hf(hf<eps)=0;
        Dh=h-hf;
        I=Dh<eps;
        Dh(I)=0;
        N=rho.*g.*Dh;
        dNdh=rho.*g;
        dNdh(I)=0;
        
    end
    
    function [N,dNdh]=NRosier()
        
        hf=rhow.*H./rho;
        hf(hf<eps)=0;
        Dh=h-hf;
        I=Dh<eps;
        Dh(I)=0;
        
        N=rho.*g.*Dh;
        dNdh=rho.*g;
        gammaw=CtrlVar.Hydrology.gammaw ;
        
        if gammaw< 1
            
            NR=  gammaw*rho.*g.*h;
            IR=NR<N ;  % when to use Rosier
            N(IR)=NR;
            
            
            dNdh(IR)=gammaw.*rho.*g;
        end
        
        dNdh(I)=0;   % h < hf
        
    end
    
    
    
    function [taubxi,taubyi,dtaubxdui,dtaubxdvi,dtaubydui,dtaubydvi,dtaubxdhi,dtaubydhi] =Weertman()
        
        taubxi=He.*beta2i.*ub; % this is the straightforward (linear) expression for basal stress
        taubyi=He.*beta2i.*vb;
        
        % Dbeta2i=(1./m-1).*(C+CtrlVar.Czero).^(-1./m).*(ub.^2+vb.^2+CtrlVar.SpeedZero^2).^((1-3*m)./(2*m));
        
        if Nouts
            
            dtaubxdui=He.*(beta2i+Dbeta2i.*ub.*ub);
            dtaubydvi=He.*(beta2i+Dbeta2i.*vb.*vb);
            
            dtaubxdvi=He.*Dbeta2i.*ub.*vb;
            dtaubydui=dtaubxdvi;
            
            dtaubxdhi=delta.*beta2i.*ub ;
            dtaubydhi=delta.*beta2i.*vb ;
            
        else
            
            dtaubxdui=[];
            dtaubxdvi=[];
            dtaubydui=[];
            dtaubydvi=[];
            dtaubxdhi=[];
            dtaubydhi=[];
            
        end
        
    end
    
    
    function [taubxi,taubyi,dtaubxdui,dtaubxdvi,dtaubydui,dtaubydvi,dtaubxdhi,dtaubydhi] =Budd()
        
        %
        % taux = G  N^(q/m) beta2 u
        %
        % dtaux/du = G N^(q/m) (beta2 + u dbeta2/du)  
        % 
        % dtaux/dh= G (q/m) N^(q/m-1) dN/dh  beta2 u  + dG/dh G  N^(q/m) beta2 u 
        %
        qm=q./m;
        Nqm=N.^(qm) ;
        
        T=He.* Nqm.*beta2i ;
        taubxi=T.*ub ;
        taubyi=T.*vb ;

        if Nouts>2
            
            dtaubxdui=He.*Nqm.*(beta2i+Dbeta2i.*ub.*ub);
            dtaubydvi=He.*Nqm.*(beta2i+Dbeta2i.*vb.*vb);
            
            dtaubxdvi=He.*Nqm.*Dbeta2i.*ub.*vb;
            dtaubydui= dtaubxdvi ;
            
            %
            % E=qm.*N.^(qm-1).*rho.*g.*(delta.*Dh+He).*beta2i ;
            %
            %  
            
            E=delta.*T+He.*qm.*N.^(qm-1).*dNdh.*beta2i ;
            
            dtaubxdhi=  E.*ub;
            dtaubydhi=  E.*vb;
            %dtaubxdhi=Nqm.*qm.*beta2i.*ub./(Dh+eps)  + qm.*delta.*N.^(qm-1).*rho.*g.*Dh.*beta2i.*ub ;
            %dtaubydhi=Nqm.*qm.*beta2i.*vb./(Dh+eps)  + qm.*delta.*N.^(qm-1).*rho.*g.*Dh.*beta2i.*vb ;
            
        else
            
            dtaubxdui=[];
            dtaubxdvi=[];
            dtaubydui=[];
            dtaubydvi=[];
            dtaubxdhi=[];
            dtaubydhi=[];
            
        end
        
    end
    
    % nested functions
    function [taubxi,taubyi,dtaubxdui,dtaubxdvi,dtaubydui,dtaubydvi,dtaubxdhi,dtaubydhi] =Coulomb()
        % Coulomb
        % tx=muk rho g (h-hf)    u/speed
        % tv=muk rho g (h-hf)    v/speed
        
        speed=sqrt(ub.*ub+vb.*vb+CtrlVar.SpeedZero^2) ;
        Tau=muk.*N ;  % muk rho g (h-hf)
        taubxi=Tau.*ub./speed;
        taubyi=Tau.*vb./speed;
        
        
        if Nouts>2
            
            temp=speed.^3 ; % (ub.*ub+vb.*vb+CtrlVar.SpeedZero^2).^{3/2) ;
            
            dtaubxdui= Tau.*( 1./speed-ub.^2./temp) ;
            dtaubydvi= Tau.*(1./speed-vb.^2./temp );
            dtaubxdvi=-Tau.*ub.*vb./temp ;
            dtaubydui=dtaubxdvi;  % just symmetry, always true,
            
            % dN/dh = rho g 
            E=muk.*dNdh./speed ;
            E(Dh<eps)=0;
            
            dtaubxdhi =  E.*ub;
            dtaubydhi =  E.*vb;
            
        else
            
            dtaubxdui=[];
            dtaubxdvi=[];
            dtaubydui=[];
            dtaubydvi=[];
            dtaubxdhi=[];
            dtaubydhi=[];
            
        end
        
    end
    
    
end