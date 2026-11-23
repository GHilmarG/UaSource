function [p,UserVar,RunInfo]=UaOptimisation(UserVar,CtrlVar,RunInfo,MUA,func,p,plb,pub,HessianSetup,InverseFileRestartWriter)

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
narginchk(8,10)
if nargin<10
    InverseFileRestartWriter=[];   % (23 Nov 2026) function handle for writing restart files during the inversion (InvertForModelParameters.m)
end
nargoutchk(3,3)


switch CtrlVar.Inverse.MinimisationMethod

    case  "-UaOptimisation-HessianBased-"


        [p,UserVar,RunInfo]=UaOptimisationHessianBased(UserVar,CtrlVar,RunInfo,MUA,func,p,plb,pub);


    case  "-UaOptimisation-HessianVectorBased-"

        % (15 November 2026) inexact Newton with Hessian-vector products, see UaOptimisationHessianVectorBased.m
        if nargin<9 || isempty(HessianSetup)
            error("UaOptimisation:NoHessianSetup","-UaOptimisation-HessianVectorBased- requires the function handle HessianSetup (see InvertForModelParameters.m).")
        end
        [p,UserVar,RunInfo]=UaOptimisationHessianVectorBased(UserVar,CtrlVar,RunInfo,MUA,func,p,plb,pub,HessianSetup,InverseFileRestartWriter);

    case  "-UaOptimisation-GradientBased-"

        if  isempty(plb)  && isempty(pub)

       
            [p,UserVar,RunInfo]=UaOptimisationGradientBased(UserVar,CtrlVar,RunInfo,MUA,func,p,plb,pub,InverseFileRestartWriter) ;

        else

            fprintf("Box contraints are enforced using bounded gradients. \n ")
            [p,UserVar,RunInfo]=UaOptimisationGradientBasedBounded(UserVar,CtrlVar,RunInfo,MUA,func,p,plb,pub,InverseFileRestartWriter) ;

        end


    otherwise

        error("CaseNotFound")


end


end