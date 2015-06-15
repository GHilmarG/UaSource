function  [Experiment,MUA,BCs,s,b,h,S,B,ub,vb,ud,vd,l,alpha,rho,rhow,g,GF,InvStartValues,Priors,Meas,BCsAdjoint,Info]=...
    GetInputsForInverseRestartRun(CtrlVar)

fprintf(CtrlVar.fidlog,' loading adjoint restart file: %s \t ',CtrlVar.NameOfAdjointRestartFiletoRead);


load(CtrlVar.NameOfAdjointRestartFiletoWrite,...
    'Experiment','CtrlVar','MUA','BCs','s','b','h','S','B','ub','vb','ud','vd','l','alpha','rho','rhow','g','GF',...
    'InvStartValues','Priors','Meas','BCsAdjoint','Info','InvFinalValues');

fprintf(CtrlVar.fidlog,' done \n ');


% Set start values to last estimates
InvStartValues=InvFinalValues;

end


