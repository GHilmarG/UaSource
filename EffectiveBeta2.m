function [taux,tauy,dtauxdu,dtauxdv,dtauydu,dtauydv,dtauxdh,dtauydh] = EffectiveBeta2(CtrlVar,He,delta,s,b,S,B,rho,rhow,ub,vb,C,m,uo,vo,Co,mo,ua,va,Ca,ma)
    
    %      He = HeavisideApprox(kH,h-hf,CtrlVar.Hh0); 
    %   delta = DiracDelta(kH,h-hf,CtrlVar.Hh0);
    %      hf = rhow*H./rho;
    %       H = S-B 
    %
    %
    % Returns basal drag and the directional derivatives of basal drag with respect to u,
    % v and h.
    %
    %
    %
    %
   
    %Cint(Cint<CtrlVar.Cmin)=CtrlVar.Cmin; % for higher order elements it is possible for Cint to be less than zero  even for C strictly positiv
    %beta2int=Cint.^(-1/m).*(sqrt(uint.*uint+vint.*vint+CtrlVar.SpeedZero^2)).^(1/m-1) ;
    
    if nargin<13
        U=u; 
        V=v;
    else
        U=u-uo;
        V=v-vo;
    end
    
    % To do: 
    %
    % Allow for different slipperiness parameters for 
    % beta2= Heint 
    
    
    beta2=He.*(C+CtrlVar.Czero).^(-1./m).*(sqrt(U.*U+V.*V+CtrlVar.SpeedZero^2)).^(1./m-1) ;
    
    taux=beta2.*u;
    tauy=beta2.*v;
    
    if nargout>1
        
        % The directional derivative is
        % D beta^2(u,v)[Delta u, Delta v]= (1/m-1) C^(-1/m) (u^2+v^2)^((1-3m)/2m)  (u \Delta u + v \Delta v)
        %
        
        %Dbeta2int=(1/m-1).*Cint.^(-1/m).*(U.^2+V.^2+CtrlVar.SpeedZero^2).^((1-3*m)/(2*m));
        Dbeta2=He.*(1./m-1).*(C+CtrlVar.Czero).^(-1./m).*(U.^2+V.^2+CtrlVar.SpeedZero^2).^((1-3*m)./(2*m));
        
        
        dtauxdu=beta2+Dbeta2.*U.*U;
        dtauydv=beta2+Dbeta2.*V.*V;
        
        dtauxdv=Dbeta2.*U.*V;
        dtauydu=dtauxdv;
        
        if any(isnan(Dbeta2)) ; save TestSave  ;  error(' NaN in Dbeta2int ' ) ; end
        
    end
    
    if nargout>5
    
       dtauxdh=delta.*(C+CtrlVar.Czero).^(-1./m).*(U.^2+V.^2+CtrlVar.SpeedZero^2).^((1./m);
       dtauydh=dtauxdh;  
       
    end
    
    
    if any(isnan(beta2)) ; save TestSave  ;  error(' NaN in beta2int ' ) ; end
    
end