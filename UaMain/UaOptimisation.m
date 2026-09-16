function [p,UserVar,RunInfo]=UaOptimisation(UserVar,CtrlVar,RunInfo,MUA,func,p,plb,pub)

%
% func is the function to me minimized
%  p is the parameter set, i.e. func(p)
%
%  Func is func evaluated as a function of step-size gamma in the direction of
%  the gradient: Func=@(gamma) func(p-gamma*dJdp);
%
%
% A very short and concise overview over some of the ideas used is found in:
%
% https://www.epfl.ch/labs/anchp/wp-content/uploads/2018/05/part5-1.pdf
%
%%
narginchk(8,8)
nargoutchk(3,3)


switch CtrlVar.Inverse.MinimisationMethod

    case  "-UaOptimization-HessianBased-"


        [p,UserVar,RunInfo]=UaOptimisationHessianBased(UserVar,CtrlVar,RunInfo,MUA,func,p,plb,pub);


    case  "-UaOptimization-GradientBased-"


       %[p,UserVar,RunInfo]=UaOptimisationGradientBased(UserVar,CtrlVar,RunInfo,MUA,func,p,plb,pub) ;
       [p,UserVar,RunInfo]=UaOptimisationGradientBasedBounded(UserVar,CtrlVar,RunInfo,MUA,func,p,plb,pub) ;


    otherwise

        error("CaseNotFound")


end


end