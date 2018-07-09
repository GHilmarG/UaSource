function [tbx,tby,tb,beta2] = CalcBasalTraction(CtrlVar,MUA,ub,vb,C,m,GF)

%
%    [tbx,tby,tb,beta2] = CalcBasalTraction(CtrlVar,MUA,ub,vb,C,m,GF)
% 
% Calculates basal traction from basal velocity using the sliding law.
% 
% Returns nodal values 
%
% Note: This can only be used to calculate basal traction when using the SSTREAM
% and the Hybrid flow approximation. This will not return correct results for
% the SSHEET approximation!
%
% Note: There is a slight inconsistency with respect to how this is done
% internally in Ua in the sense that the floating mask is here evalutated at
% nodes, whereas internally this is done at integration pointa. 
%
%
%

narginchk(7,7)

if CtrlVar.CisElementBased
    % project onto nodes
    [M,ElePerNode] = Ele2Nodes(MUA.connectivity,MUA.Nnodes); 
    
    C=M*C;
    m=M*m;
    
end


beta2=(C+CtrlVar.Czero).^(-1./m).*(sqrt(ub.*ub+vb.*vb+CtrlVar.SpeedZero^2)).^(1./m-1) ;

tbx=GF.node.*beta2.*ub; 
tby=GF.node.*beta2.*vb;
tb=sqrt(tbx.^2+tby.^2);


end

