function [dJdpModified,RunInfo]=NextGradient(dJdp,dJdpLast,dJdpModified,CtrlVar,RunInfo)


switch lower(CtrlVar.GradientUpgradeMethod)
    
    case 'conjgrad'
        
        [dJdpModified,RunInfo]=NewConjugatedGrad(-dJdp,-dJdpLast,-dJdpModified,CtrlVar,RunInfo);
        
    otherwise
   
        dJdpModified=dJdp;
end



end