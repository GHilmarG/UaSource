
function Hessian = BuildInversionHessian(CtrlVar,MUA,F,BCs,l,Priors,Meas,BCsAdjoint,Psi_x,Psi_y)


% Typically, in any large-scale inversions, the Hessian is never calculated here. Rather, an Hessian approximation is build
% up using a BFSG update.
%
% However, for some smaller sized problems the Hessian, or some approximations thereof, can be calculated and returned.
% 
% There are currently three options
%
% 1) Calculate the Hessian using the direct-adjoint method. This gives a full Hessian and a very exact estimate. This is as
% good as it gets, but requires large memory and gives a full Hessian!  This will only work for small to medium sized
% problems with a few tens of thousand nodes.
%
% 2) Calculate the Hessian using finite-differences as done by the MATLAB toolbox. This is done in a highly optimized way and
% works surprisingly good. For some reason the MATLAB toolbox wants a zero sparse matrix returned. And this is done here, but
% the Hessian approximation is calculated by the toolbox function.
%
% 3) Some ad-hoc approximations of the Hessian. This could be as simple as using the Metric Matrix or something similar,
% Sobolev smoothing, etc. None of these approaches work particularly well, and were here implemented basically out of
% a mixture of curiosity and unjustified optimism.
%


switch CtrlVar.Inverse.Hessian

    case "-DirectAdjoint-"

        Hessian = CalcDirectAdjointHessian(CtrlVar,MUA,F,BCs,l,Priors,Meas,BCsAdjoint,Psi_x,Psi_y) ;

    case "-Jpp-"

        Hessian=Jpp(CtrlVar,MUA,Meas);

    case "-FiniteDifferences-"


        % Not sure I have fully correctly understood this, but it appears that the Matlab toolbox functions might still need a
        % Hessian argument returned, even if the Hessian is calculated using finite differences. So here I am returning an empty
        % sparse matrix of the right size.
        [isA,isB,isC] = isABC(CtrlVar);
        np=(isA+isB+isC)*MUA.Nnodes;
        Hessian=sparse(np,np);

    otherwise

        warning("BuildInversionHessian:CaseNotFound","While building inverse Hessian, case was not found. Returning empty matrix")

        Hessian=[];
end

end