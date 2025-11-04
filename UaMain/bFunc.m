function [J,f,H]=bFunc(b,CtrlVar,s,B,S,rho,rhow)


hf=rhow*(S-B)./rho ;
h=s-b; 
G = HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);
dGdb=-DiracDelta(CtrlVar.kH,h-hf,CtrlVar.Hh0) ;


f=b - G.*B - (1-G).*(rho.*s-rhow.*S)./(rho-rhow) ;
J=sum(f.^2)/2 ;


if nargout> 2
    elements = 1 - dGdb.* (B -  (rho.*s-rhow.*S)./(rho-rhow)) ;
    I=1:numel(b) ; 
    H=sparse(I,I,elements) ;

end




end