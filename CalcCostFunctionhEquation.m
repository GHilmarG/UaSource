function [r,UserVar,RunInfo,rForce,rWork,D2]=CalcCostFunctionhEquation(UserVar,RunInfo,CtrlVar,MUA,gamma,F1,F0,L,Lrhs,l,dh,dl)



narginchk(12,12)
nargoutchk(1,6)


F1.h=F1.h+gamma*dh;
l=l+gamma*dl;


% Make sure to update all other fields that depend on h
% Here the only such field is potentially the mass balance term

 CtrlVar.ResetThicknessToMinThickness=0;
 [F1.b,F1.s]=Calc_bs_From_hBS(CtrlVar,MUA,F1.h,F1.S,F1.B,F1.rho,F1.rhow);
 [UserVar,F1]=GetMassBalance(UserVar,CtrlVar,MUA,F1); % actually this call only needed if mass-balance depends on h

% Only here evaluating the righ-hand side of the equation, is the J(x0+ gamma dx)
[UserVar,R]=MassContinuityEquationAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1) ;


if ~isempty(L)
    
    frhs=-R-L'*l;       
    grhs=Lrhs-L*F1.h;  
    
else
    frhs=-R;
    grhs=[];
end

frhs=frhs/MUA.Area;  % It's OK to do this only here, because I scale frhs and grhs equally.
grhs=grhs/MUA.Area;


D2=[frhs;grhs]'*[dh;dl];
rWork=full(D2^2);

rForce=full([frhs;grhs]'*[frhs;grhs]);



switch CtrlVar.hMinimisationQuantity
    case "Force Residuals"
        r=rForce;
    case "Work Residuals"
        r=rWork;
end


end

