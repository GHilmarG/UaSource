function [gradJ1,RunInfo]=NextGradient(dJdp,dJdpLast,gradJ0,CtrlVar,RunInfo)


switch lower(CtrlVar.GradientUpgradeMethod)
    
    case 'conjgrad'
        
        [gradJ1,RunInfo]=NewConjugatedGradMetric(dJdescent,dJdescentlast,gradJ0,G,CtrlVar,RunInfo);
        [dJdpModified,RunInfo]=NewConjugatedGrad(-dJdp,-dJdpLast,-dJdpModified,CtrlVar,RunInfo);
        
    otherwise
   
        dJdpModified=dJdp;
end



end