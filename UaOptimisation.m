function [p,UserVar,RunInfo]=UaOptimisation(UserVar,CtrlVar,RunInfo,MUA,func,p,plb,pub)

%
% func is the function to me minimized
%  p is the paramter set, i.e. func(p)
%
%  Func is func evaluated as a function of stepsize gamma in the direction of
%  the gradient: Func=@(gamma) func(p-gamma*dJdp);
%

narginchk(8,8)
nargoutchk(3,3)



switch CtrlVar.Inverse.MinimisationMethod
    
    
    case {"UaOptimization","UaOptimization-Gradient"}
        
        [p,UserVar,RunInfo]=UaOptimisationGradientBased(UserVar,CrlVar,RunInfo,MUA,func,p,plb,pub) ;
        
        
    case "UaOptimization-Hessian"
        
        [p,UserVar,RunInfo]=UaOptimisationHessianBased(UserVar,CtrlVar,RunInfo,MUA,func,p,plb,pub) ;
        
    otherwise
        
        
        error('Case not found')
        
        
end


end