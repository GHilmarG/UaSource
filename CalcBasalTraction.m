function [tbx,tby,tb,beta2] = CalcBasalTraction(CtrlVar,MUA,ub,vb,C,m,GF)

%[tbx,tby,tb,beta2] = CalcBasalTraction(CtrlVar,MUA,ub,vb,C,m,GF)
% calculates basal traction from basal velocity using the sliding law
% returns nodal values 
% t_{bx}= C^{-1/m} |v|^{1/m-1} u
% t_{bv}= C^{-1/m} |v|^{1/m-1} v
% |v|=C \tau_b^m

narginchk(7,7)

if CtrlVar.CisElementBased
    % project onto nodes
    [M,ElePerNode] = Ele2Nodes(MUA.connectivity,MUA.Nnodes); C=M*C;
    
end

% I calculate: 
% beta2int=(Cint+CtrlVar.Czero).^(-1/m).*(sqrt(uint.*uint+vint.*vint+CtrlVar.SpeedZero^2)).^(1/m-1) ;
% therefore both speed and C must be modified accordingly 
% taub=beta2 ub
% In the limiting case C=0 and ub=0 I get
% beta2=Czero^(-1/m) ubzero^(2/m-2)
% 
%
%
C=C+CtrlVar.Czero;
speed=sqrt(real(ub.*ub+vb.*vb+CtrlVar.SpeedZero^2));

tb=real((speed./C).^(1./m));
tbx=real(C.^(-1./m)).*real(speed.^(1./m-1)).*ub;
tby=real(C.^(-1./m)).*real(speed.^(1./m-1)).*vb;

tb=tb.*GF.node ; tbx=tbx.*GF.node ; tby=tby.*GF.node ;

beta2=C.^(-1./m).*speed.^(1./m-1);

end

