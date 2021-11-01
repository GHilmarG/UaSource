function [r,UserVar,RunInfo,rForce,rWork,D2]=CalcCostFunctionhEquation(UserVar,RunInfo,CtrlVar,MUA,gamma,F1,F0,L,Lrhs,l,dh,dl)



narginchk(12,12)
nargoutchk(1,6)


F1.h=F1.h+gamma*dh;
l=l+gamma*dl;


[UserVar,R]=MassContinuityEquationAssembly(UserVar,CtrlVar,MUA,F0.h,F0.rho,F0.ub,F0.vb,F0.as,F0.ab,F1.h,F1.ub,F1.vb,F1.as,F1.ab,F1.dasdh,F1.dabdh);

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

