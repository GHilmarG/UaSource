function [r,UserVar,RunInfo,rForce,rWork,D2]=CalcCostFunctionLevelSetEquation(UserVar,RunInfo,CtrlVar,MUA,gamma,F1,F0,L,Lrhs,l,dLSF,dl,BCs,R)

%
% (25 Oct 2026) Optional input R: the residual at gamma=0, if already assembled (in LevelSetEquationNewtonRaphson.m). It is only
% used for gamma=0, where it is the residual at F1, and the residual is then not assembled again.
%

narginchk(13,14)
nargoutchk(1,6)

%
%
%
%

MLC=BCs2MLC(CtrlVar,MUA,BCs);
L=MLC.LSFL ; Lrhs=MLC.LSFRhs ;

if isempty(l) || (numel(l)~=numel(Lrhs))
    l=Lrhs*0 ;
end

if isempty(dl) || (numel(dl)~=numel(Lrhs))
    dl=Lrhs*0 ;
end

if isempty(dLSF) || (numel(dLSF)~=numel(F1.LSF))
    dLSF=F1.LSF*0;
end


F1.LSF=F1.LSF+gamma*dLSF;
l=l+gamma*dl;

if nargin<14 || isempty(R) || gamma~=0
    [UserVar,R]=LevelSetEquationAssemblyNR2(UserVar,CtrlVar,MUA,F0,F1) ;
end
 %[UserVar,R]=LevelSetEquationAssemblyNR2(UserVar,CtrlVar,MUA,F0.LSF,F0.c,F0.ub,F0.vb,F1.LSF,F1.c,F1.ub,F1.vb,F0.LSFqx,F0.LSFqy,F1.LSFqx,F1.LSFqy);


if ~isempty(L)

    frhs=-R-L'*l;        % Units: Area  [\varphi] / time
    grhs=Lrhs-L*F1.LSF;  % Units: [\varphi]=distance, but this should always be zero anyhow if initial point is feasible


else
    frhs=-R;
    grhs=[];
end

frhs=frhs/MUA.Area;  % It's OK to do this only here, because I scale frhs and grhs equally.
grhs=grhs/MUA.Area;


D2=[frhs;grhs]'*[dLSF;dl];
rWork=full(D2^2);

rForce=full([frhs;grhs]'*[frhs;grhs]);



switch CtrlVar.LSFMinimisationQuantity
    case "Force Residuals"
        r=rForce;
    case "Work Residuals"
        r=rWork;
end


end

