function [N,dNdh]=EffectivePressure(CtrlVar,peffModel,s,b,S,B,h,rho,rhow,g)

narginchk(9,9)

switch   peffModel

    case "N0"

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

    case "N1"


    otherwise


        error(" peffModel?")

end





end
