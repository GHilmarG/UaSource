function  [MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,alpha,rho,rhow,g,GF,InvStartValues,Priors,Meas,BCsAdjoint,Info]=...
    GetInputsForInverseRestartRun(CtrlVar)

fprintf(CtrlVar.fidlog,' loading adjoint restart file: %s \t ',CtrlVar.NameOfAdjointRestartFiletoRead);


load(CtrlVar.NameOfAdjointRestartFiletoRead,...
    'MUA','BCs','s','b','h','S','B','ub','vb','ud','vd','l','alpha','rho','rhow','g','GF',...
    'InvStartValues','Priors','Meas','BCsAdjoint','Info','InvFinalValues');

fprintf(CtrlVar.fidlog,' done \n ');


% Set start values to last estimates
InvStartValues=InvFinalValues;


[InvStartValues.AGlen,InvStartValues.n]=TestAGlenInputValues(CtrlVar,MUA,InvStartValues.AGlen,InvStartValues.n);
[Priors.AGlen,Priors.n]=TestAGlenInputValues(CtrlVar,MUA,Priors.AGlen,Priors.n);

[InvStartValues.C,InvStartValues.m]=TestSlipperinessInputValues(CtrlVar,MUA,InvStartValues.C,InvStartValues.m);
[Priors.C,Priors.m]=TestSlipperinessInputValues(CtrlVar,MUA,Priors.C,Priors.m);

[Priors.rho,Priors.rhow]=TestDensityInputValues(CtrlVar,MUA,Priors.rho,Priors.rhow);




end


