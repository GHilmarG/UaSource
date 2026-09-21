




function FppFiniteDifferences(CtrlVar,MUA,F,BCs,l,BCsAdjoint,Psi_x,Psi_y)


%% log c test

CtrlVar.log10Derivatives=true;
%F=PsiTimesddFuvdCdC(CtrlVar,MUA,F,Psi_x,Psi_y); 
HFCC=FCC(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);

C0=F.C; 
logC0=log10(C0);

iColumn=7;  % just do the finite-difference comparison for this column of the Hessian contribution Fpp.


CtrlVar.Inverse.TestAdjoint.FiniteDifferenceRelStepSize=0.0001;

% comparison 

F.C=C0; 
% perturbation 

if CtrlVar.log10Derivatives
    % perturb in log(C) space
  
    dlogC=1e-5;
   
else
    % perturb in linear space
    deltaStep=CtrlVar.Inverse.TestAdjoint.FiniteDifferenceRelStepSize*abs(C0)+CtrlVar.Inverse.TestAdjoint.FiniteDifferenceStepSize;
    dC=deltaStep(iColumn);
end

% comparison 

if CtrlVar.log10Derivatives
    F.C(iColumn)=10^(logC0(iColumn)-dlogC);
else
    F.C(iColumn)=F.C(iColumn)-dC;
end

FpMinus=dIdCq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y); 

F.C=C0; 

if CtrlVar.log10Derivatives
   F.C(iColumn)=10^(logC0(iColumn)+dlogC);
else
    F.C(iColumn)=F.C(iColumn)+dC;
end

FpPlus=dIdCq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y); 

if CtrlVar.log10Derivatives
    HFCC_col_FD=(FpPlus-FpMinus)/(2*dlogC);  % Finite difference estimate of the iColumn of the Hessian

else
    HFCC_col_FD=(FpPlus-FpMinus)/(2*dC);
end



norm(HFCC(:,iColumn)-HFCC_col_FD)/norm(HFCC_col_FD)   


FTest=FindOrCreateFigure("Fpp Test") ; plot(HFCC(:,iColumn),HFCC_col_FD,"or") ; axis equal ; 
hold on ; 
plot([min(HFCC_col_FD) max(HFCC_col_FD)],[min(HFCC_col_FD) max(HFCC_col_FD)],"--k")  

% ax=gca ; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin'; axis on ; axis equal tight ; box off 

xlabel("Direct-Adjoint",Interpreter="latex")  ;
ylabel("Finite difference",Interpreter="latex")
title("$\mathcal{F}^{pp}_{CC,lm} = \langle \Psi_x,\delta^2_{CC}\mathcal{F}_x[\phi_l,\phi_m] \rangle  + \langle \Psi_y , \delta^2_{CC}\mathcal{F}_y[\phi_l,\phi_m] \rangle$",Interpreter="latex")
subtitle("Comparison is here for one random column",Interpreter="latex")


discrepancy = HFCC_col_FD(iColumn) - HFCC(iColumn,iColumn);         % true minus your formula's diagonal

ln10 = log(10);



b_C=dIdCq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);
correction_used = b_C(:)*ln10;               % the term you added
ratio = discrepancy ./ correction_used;


fprintf("                                                  AD                        FD \n")
n=10 ; [(1:n)' full(HFCC(1:n,iColumn)) HFCC_col_FD(1:n) full(HFCC(1:n,iColumn))-HFCC_col_FD(1:n)]  
fprintf("correction used at (%i,%i) is %g \n",iColumn,iColumn,correction_used(iColumn))

ratio(1:n)

end 
