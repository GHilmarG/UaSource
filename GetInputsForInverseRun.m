function [UserVar,InvStartValues,Priors,Meas,BCsAdjoint]=GetInputsForInverseRun(UserVar,CtrlVar,MUA,BCs,F,l,GF)

%[UserVar,InvStartValues,Priors,Meas,BCsAdjoint]=GetInputsForInverseRun(UserVar,CtrlVar,MUA,BCs,time,AGlen,C,n,m,s,b,S,B,rho,rhow,GF,g,alpha,ub,vb,ud,vd,l)



BCsAdjoint=BoundaryConditions;
Meas=Measurements;
Priors=PriorProbabilityDistribution;
InvStartValues=InversionValues;

[UserVar,InvStartValues,Priors,Meas,BCsAdjoint]=DefineInputsForInverseRun(UserVar,CtrlVar,MUA,BCs,InvStartValues,Priors,Meas,BCsAdjoint,CtrlVar.time,F.AGlen,F.C,F.n,F.m,F.s,F.b,F.S,F.B,F.rho,F.rhow,GF,F.g,F.alpha,F.ub,F.vb,F.ud,F.vd,l);

BCsAdjoint=CreatePlausibleBCsForAdjointProblem(BCs,BCsAdjoint);


%% Test inputs
Meas=TestMeas(CtrlVar,MUA,Meas);

[InvStartValues.AGlen,InvStartValues.n]=TestAGlenInputValues(CtrlVar,MUA,InvStartValues.AGlen,InvStartValues.n);
[Priors.AGlen,Priors.n]=TestAGlenInputValues(CtrlVar,MUA,Priors.AGlen,Priors.n);

[InvStartValues.C,InvStartValues.m]=TestSlipperinessInputValues(CtrlVar,MUA,InvStartValues.C,InvStartValues.m);
[Priors.C,Priors.m]=TestSlipperinessInputValues(CtrlVar,MUA,Priors.C,Priors.m);

[Priors.rho,Priors.rhow]=TestDensityInputValues(CtrlVar,MUA,Priors.rho,Priors.rhow);

isCorrectDimentions=DoPriorsHaveCorrectDimentions(CtrlVar,MUA,Priors);
if ~ isCorrectDimentions
    fprintf(' Priors do not have right dimentions at restart. \n')
    fprintf(' Modify DefineInputsForInverseRun to ensure that dimentions are correct.\n')
    error('Ua:GetInputForInverseRun:incorrectdimentisons','incorrect dimentions')
    
end

if isempty(InvStartValues.AGlen) ; save TestSave ; error('GetInputsForInverseRun:empty','InvStartValues.AGlen is empty') ; end
if isempty(InvStartValues.C) ; save TestSave ; error('GetInputsForInverseRun:empty','InvStartValues.C is empty') ; end
if isempty(InvStartValues.n) ; save TestSave ; error('GetInputsForInverseRun:empty','InvStartValues.n is empty') ; end
if isempty(InvStartValues.m) ; save TestSave ; error('GetInputsForInverseRun:empty','InvStartValues.m is empty') ; end


end

























