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

    % rounding and internal element interpolation can make these negative.
    % If not made positive, some variables can be complex numbers on
    % return.
    delta(delta<0)=0;
    He(He<0)=0 ;


    %% Basal drag term : ice
    % this drag term is zero if the velocities are zero.
    
    C0=CtrlVar.Czero;
    u0=CtrlVar.SpeedZero;
    
    speed=(sqrt(ub.*ub+vb.*vb+u0^2)); 
    Um=speed.^(1./m-1) ;
    beta2i=(C+C0).^(-1./m).*Um ; %   (sqrt(ub.*ub+vb.*vb+CtrlVar.SpeedZero^2)).^(1./m-1) ;

    
    % Dbeta2i is zero for m=1.
    Dbeta2i=(1./m-1).*(C+C0).^(-1./m).*(ub.^2+vb.^2+u0^2).^((1-3*m)./(2*m));
    
    
    if ~isfield(CtrlVar,"SlidingLaw")
        CtrlVar.SlidingLaw="Weertman";
    end
    
    if ~isfield(CtrlVar.Inverse,'dFuvdClambda')
        CtrlVar.Inverse.dFuvdClambda=false;
    end
    
    if CtrlVar.Inverse.dFuvdClambda
        
        % dF/dC
        %  Just take the derivative of the tau term with respect to C. And don't forget
        %  the minus in front of the tau term in the momemtum equation.
        %
        %  -Taux = |Tau|  u/U, 
        %  -Tauy = |Tau|  v/U, 
        %   where U is the speed.
        %
        % I include the u and v in the adjoint calculation itself, so I just need the
        % derivative:
        %
        %   d (-|Tau|/U) / dC 
        
       
        
        switch CtrlVar.SlidingLaw
            case {"W","Weertman"}
                
                % tau = He * (C+CtrlVar.Czero).^(-1./m)   * U
                % 
                % (U^(1/m - 1)*He(h - hf))/(m*(C + C0)^(1/m + 1))
                
                dFuvdC =  He.*    (1./m).*(C+C0).^(-1./m-1)   .*Um;
                
            case {"B","Budd","W-N0"}
    
                %
                %             dFuvdC= Nqm.*(1./m).*(C+C0).^(-1./m-1)  .*U;
                %  just take the derivative with respect to C
                %   tau = He.*Nqm.* beta2i
                %       = He.*Nqm.* (C+C0.^(-1./m).*U ; %   with U=(sqrt(ub.*ub+vb.*vb+CtrlVar.SpeedZero^2)).^(1./m-1) ;
                %
                
                hf=rhow.*H./rho;
                hf(hf<eps)=0;
                Dh=h-hf; Dh(Dh<eps)=0;
                N=rho.*g.*Dh;
                qm=q./m;
                Nqm=N.^(qm) ;
                
                
                dFuvdC= He.*Nqm.*(1./m).*(C+C0).^(-1./m-1)  .*Um;
                
            case {"rpCW-N0","Cornford"}
                
                U=speed;
                N=N0(CtrlVar,h,H,rho,rhow,g) ;
                dFuvdC=(U.^(1.0./m-1.0).*muk.^(m+1.0).*N.^(m+1.0).*He.*(muk.^m.*N.^m+U.*(He.*(C+C0).^(-1.0./m)).^m).^(-1.0./m-1.0).*(C+C0).^(-1.0./m-1.0))./m;
                
            case {"rCW-N0","Umbi"} % reciprocal Coulumb-Weertman with zeroth-order hydrology
                
                U=speed;
                N=N0(CtrlVar,h,H,rho,rhow,g) ;
                dFuvdC=(U.^(1.0./m).*muk.^2.*N.^2.*He.*1.0./(U.^(1.0./m).*He+muk.*N.*(C+C0).^(1.0./m)).^2.*(C+C0).^(1.0./m-1.0))./(U.*m) ;
                
            case {"Tsai","minCW-N0"}
                
                dFuvdC =  He.*    (1./m).*(C+C0).^(-1./m-1)   .*Um;

                N=N0(CtrlVar,h,H,rho,rhow,g);
                [taubxi,taubyi] = Weertman(CtrlVar,He,delta,ub,vb,beta2i,Dbeta2i) ;
                [taubxiC,taubyiC] = Coulomb(CtrlVar,muk,ub,vb,N) ;
                
                TauWeertman2=taubxi.*taubxi+taubyi.*taubyi ;
                TauCoulomb2=taubxiC.*taubxiC+taubyiC.*taubyiC ;
                isCoulomb=TauCoulomb2 <  TauWeertman2 ;
                
                dFuvdC(isCoulomb)=0 ; 
                
                
            case {"Coulomb","-C-"}
                
                
                fprintf("Inversion using Coulomb sliding law not implemented. \n")
                error("BasalDrag:InvalidCase","Inversion using Tsai sliding law not implemented.")
                
            otherwise
                
                error("BasalDrag:InvalidCase","Unknown case.")
                
        end
        taubx=dFuvdC;
        return
    end
    
    
    switch CtrlVar.SlidingLaw
        
        case {"Weertman","W"}
            
          [taubxi,taubyi,dtaubxdui,dtaubxdvi,dtaubydui,dtaubydvi,dtaubxdhi,dtaubydhi] = Weertman(CtrlVar,He,delta,ub,vb,beta2i,Dbeta2i) ;
            
          % Weertman2 is done using symbolic toolbox, just did this out of curiosity, the
          % results are exactly the same.
          % [taubxi,taubyi,dtaubxdui,dtaubydvi,dtaubxdvi,dtaubydui,dtaubxdhi,dtaubydhi] = Weertman2(C,CtrlVar.Czero,He,delta,m,ub,vb,CtrlVar.SpeedZero) ;

            
        case {"Budd","W-N0"}
            
            [N,dNdh]=N0(CtrlVar,h,H,rho,rhow,g); 
            [taubxi,taubyi,dtaubxdui,dtaubxdvi,dtaubydui,dtaubydvi,dtaubxdhi,dtaubydhi] = Budd(CtrlVar,He,delta,ub,vb,N,dNdh,beta2i,Dbeta2i,q,m) ;
            
            
        case {"Coulomb","C"}
            

            [N,dNdh]=N0(CtrlVar,h,H,rho,rhow,g); 
            
            hf=rhow.*H./rho; hf(hf<eps)=0;
            Dh=h-hf; Dh(Dh<eps)=0;
            [taubxi,taubyi,dtaubxdui,dtaubxdvi,dtaubydui,dtaubydvi,dtaubxdhi,dtaubydhi] = Coulomb(CtrlVar,muk,ub,vb,N,dNdh,Dh) ;
            
            
        case {"Tsai","minCW-N0"}
                 
            hf=rhow.*H./rho; hf(hf<eps)=0;
            Dh=h-hf; Dh(Dh<eps)=0;
            
       
            [N,dNdh]=N0(CtrlVar,h,H,rho,rhow,g); 
            [taubxi,taubyi,dtaubxdui,dtaubxdvi,dtaubydui,dtaubydvi,dtaubxdhi,dtaubydhi] = Weertman(CtrlVar,He,delta,ub,vb,beta2i,Dbeta2i) ;
            [taubxiC,taubyiC,dtaubxduiC,dtaubxdviC,dtaubyduiC,dtaubydviC,dtaubxdhiC,dtaubydhiC] = Coulomb(CtrlVar,muk,ub,vb,N,dNdh,Dh) ;
            
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
            
            
        case {"rpCW-N0","Cornford"}
            
            [N,dNdh]=N0(CtrlVar,h,H,rho,rhow,g) ;
            [taubxi,taubyi,dtaubxdui,dtaubydvi,dtaubxdvi,dtaubydui,dtaubxdhi,dtaubydhi] =  rpCWN0(C,CtrlVar.Czero,N,dNdh,He,delta,m,muk,ub,vb,CtrlVar.SpeedZero) ;
            
        case {"rCW-N0","Umbi"} % reciprocal Coulumb-Weertman with zeroth-order hydrology
            
            [N,dNdh]=N0(CtrlVar,h,H,rho,rhow,g) ;
            [taubxi,taubyi,dtaubxdui,dtaubydvi,dtaubxdvi,dtaubydui,dtaubxdhi,dtaubydhi] =  rCWN0(C,CtrlVar.Czero,N,dNdh,He,delta,m,muk,ub,vb,CtrlVar.SpeedZero) ;
            
            
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
    
        
    
end

%% local functions
function [taubxi,taubyi,dtaubxdui,dtaubxdvi,dtaubydui,dtaubydvi,dtaubxdhi,dtaubydhi] = Coulomb(CtrlVar,muk,ub,vb,N,dNdh,Dh)
    % Coulomb
    % tx=muk rho g (h-hf)    u/speed
    % tv=muk rho g (h-hf)    v/speed
    
    Nouts=nargout ;
    
    
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



function [taubxi,taubyi,dtaubxdui,dtaubxdvi,dtaubydui,dtaubydvi,dtaubxdhi,dtaubydhi] = Weertman(CtrlVar,He,delta,ub,vb,beta2i,Dbeta2i)
    
    Nouts=nargout ;
    
    taubxi=He.*beta2i.*ub; 
    taubyi=He.*beta2i.*vb;
    
    % Dbeta2i=(1./m-1).*(C+CtrlVar.Czero).^(-1./m).*(ub.^2+vb.^2+CtrlVar.SpeedZero^2).^((1-3*m)./(2*m));
    
    
    
    if Nouts > 2
        
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




function [taubxi,taubyi,dtaubxdui,dtaubxdvi,dtaubydui,dtaubydvi,dtaubxdhi,dtaubydhi] = Budd(CtrlVar,He,delta,ub,vb,N,dNdh,beta2i,Dbeta2i,q,m)
    
    % Weertman(CtrlVar,He,delta,ub,vb,beta2i,Dbeta2i)
    % taux = G  N^(q/m) beta2 u
    %
    % dtaux/du = G N^(q/m) (beta2 + u dbeta2/du)
    %
    % dtaux/dh= G (q/m) N^(q/m-1) dN/dh  beta2 u  + dG/dh G  N^(q/m) beta2 u
    %
    
    Nouts=nargout ;
    
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

function [N,dNdh]=N0(CtrlVar,h,H,rho,rhow,g)
    
    narginchk(6,6)
    
    hf=rhow.*H./rho;
    hf(hf<eps)=0;  % positive floation thickness 
    Dh=h-hf;
    I=Dh<eps;
    Dh(I)=0;
    
    N=rho.*g.*Dh+ CtrlVar.Nzero ; 
    dNdh=rho.*g;
    dNdh(I)=0;
    
    % Testing
   %  N=rho.*g.*10;
   %  dNdh=dNdh*0; 
    
    
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























