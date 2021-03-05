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


if contains(CtrlVar.Inverse.MinimisationMethod,"Hessian")
    
    [p,UserVar,RunInfo]=UaOptimisationHessianBased(UserVar,CtrlVar,RunInfo,MUA,func,p,plb,pub) ;

    
else
    
    [p,UserVar,RunInfo]=UaOptimisationGradientBased(UserVar,CtrlVar,RunInfo,MUA,func,p,plb,pub) ;
    
    
end


end