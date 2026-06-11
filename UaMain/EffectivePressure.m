function [N,dNdh]=EffectivePressure(CtrlVar,peffModel,s,b,S,B,h,rho,rhow,g)

narginchk(10,10)

switch   peffModel

    case "N0"

        H=S-B ; 
        hf=rhow.*H./rho;
        hf(hf<eps)=0;  % positive floation thickness
        Dh=h-hf;
        I=Dh<eps;
        Dh(I)=0;   % h above flotation
        

        N=rho.*g.*Dh+ CtrlVar.Nzero ;
        dNdh=rho.*g;
        dNdh(I)=0;

        % Testing
        %  N=rho.*g.*10;
        %  dNdh=dNdh*0;

    case "N1"

       d=S-b;
       d(d<eps)=0 ; 

       N=g*(rho.*h-rhow.*d) ; 

       dNdh=g*rho ; 

    otherwise


        error(" peffModel?")

end





end
