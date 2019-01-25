function [UserVar,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo]=GetInputsForInverseRun(UserVar,CtrlVar,MUA,BCs,F,l,RunInfo)

%[UserVar,InvStartValues,Priors,Meas,BCsAdjoint]=GetInputsForInverseRun(UserVar,CtrlVar,MUA,BCs,time,AGlen,C,n,m,s,b,S,B,rho,rhow,GF,g,alpha,ub,vb,ud,vd,l)

narginchk(7,7) 
nargoutchk(6,6)


BCsAdjoint=BoundaryConditions;
Meas=Measurements;
Priors=PriorProbabilityDistribution;
InvStartValues=InversionValues;

[UserVar,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo]=DefineInputsForInverseRun(UserVar,CtrlVar,MUA,BCs,F,l,F.GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);


BCsAdjoint=CreatePlausibleBCsForAdjointProblem(BCs,BCsAdjoint);


%% Test inputs
Meas=TestMeas(CtrlVar,MUA,Meas);

[InvStartValues.AGlen,InvStartValues.n]=TestAGlenInputValues(CtrlVar,MUA,InvStartValues.AGlen,InvStartValues.n);
[Priors.AGlen,Priors.n]=TestAGlenInputValues(CtrlVar,MUA,Priors.AGlen,Priors.n);

[InvStartValues.C,InvStartValues.m]=TestSlipperinessInputValues(CtrlVar,MUA,InvStartValues.C,InvStartValues.m);
[Priors.C,Priors.m]=TestSlipperinessInputValues(CtrlVar,MUA,Priors.C,Priors.m);

[Priors.rho,Priors.rhow]=TestDensityInputValues(CtrlVar,MUA,Priors.rho,Priors.rhow);

isCorrectDimensions=DoPriorsHaveCorrectDimensions(CtrlVar,MUA,Priors);
if ~ isCorrectDimensions
    fprintf(' Priors do not have right dimensions at restart. \n')
    fprintf(' Modify DefineInputsForInverseRun to ensure that dimensions are correct.\n')
    error('Ua:GetInputForInverseRun:incorrectdimentisons','incorrect dimensions')
    
end

if isempty(InvStartValues.AGlen) ; save TestSave ; error('GetInputsForInverseRun:empty','InvStartValues.AGlen is empty') ; end
if isempty(InvStartValues.C) ; save TestSave ; error('GetInputsForInverseRun:empty','InvStartValues.C is empty') ; end
if isempty(InvStartValues.n) ; save TestSave ; error('GetInputsForInverseRun:empty','InvStartValues.n is empty') ; end
if isempty(InvStartValues.m) ; save TestSave ; error('GetInputsForInverseRun:empty','InvStartValues.m is empty') ; end

if isempty(Priors.AGlenmax) 
    Priors.AGlenmax=CtrlVar.AGlenmax;
end

if isempty(Priors.AGlenmin)
    Priors.AGlenmin=CtrlVar.AGlenmin;
end

if isempty(Priors.Cmax)
    Priors.Cmax=CtrlVar.Cmax;
end


if isempty(Priors.Cmin)
    Priors.Cmin=CtrlVar.Cmin;
end

if isempty(Priors.Bmax)
    Priors.Bmax=Meas.s-CtrlVar.ThickMin;
end


if isempty(Priors.Bmin)
    Priors.Bmin=-1e10;
end

end

























