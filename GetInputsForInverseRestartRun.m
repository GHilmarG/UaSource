function [UserVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo]=...
    GetInputsForInverseRestartRun(UserVar,CtrlVar,RunInfo)

narginchk(3,3) 
nargoutchk(10,10)


fprintf(CtrlVar.fidlog,' Inverse run: loading restart file: %s \t ',CtrlVar.Inverse.NameOfRestartInputFile);

load(CtrlVar.Inverse.NameOfRestartInputFile,...
    'CtrlVarInRestartFile','UserVarInRestartFile','MUA','BCs','F','GF','l','RunInfo',...
    'InvStartValues','Priors','Meas','BCsAdjoint','InvFinalValues');

fprintf(CtrlVar.fidlog,' done \n ');

F.GF=GF;
F.time=CtrlVar.time ; 

LastRunInfo=RunInfo;

RunInfo=UaRunInfo;  % Hilmar check
RunInfo.Inverse=LastRunInfo.Inverse;
RunInfo.File.Name=CtrlVar.Experiment+"-RunInfo.txt";
RunInfo.File.fid = fopen(RunInfo.File.Name,'a');
% RunInfo=Validate(RunInfo);  % hilmar check
% Set start values to last estimates
InvStartValues=InvFinalValues;


[InvStartValues.AGlen,InvStartValues.n]=TestAGlenInputValues(CtrlVar,MUA,InvStartValues.AGlen,InvStartValues.n);
[InvStartValues.C,InvStartValues.m,InvStartValues.q,InvStartValues.muk]=TestSlipperinessInputValues(CtrlVar,MUA,InvStartValues.C,InvStartValues.m,InvStartValues.q,InvStartValues.muk);

%[Priors.AGlen,Priors.n]=TestAGlenInputValues(CtrlVar,MUA,Priors.AGlen,Priors.n);
%[Priors.C,Priors.m]=TestSlipperinessInputValues(CtrlVar,MUA,Priors.C,Priors.m);
%[Priors.rho,Priors.rhow]=TestDensityInputValues(CtrlVar,MUA,Priors.rho,Priors.rhow);


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




isCorrectDimensions=DoPriorsHaveCorrectDimensions(CtrlVar,MUA,Priors);


if ~isCorrectDimensions
    
    [UserVar,~,Priors,~,BCsAdjoint,RunInfo]=DefineInputsForInverseRun(UserVar,CtrlVar,MUA,BCs,F,l,F.GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);
    %[UserVar,~,Priors,~,~]=DefineInputsForInverseRun(UserVar,CtrlVar,MUA,BCs,InvStartValues,Priors,Meas,BCsAdjoint,CtrlVar.time,AGlen,C,n,m,s,b,S,B,rho,rhow,GF);
    
end

isCorrectDimensions=DoPriorsHaveCorrectDimensions(CtrlVar,MUA,Priors);
if ~ isCorrectDimensions
    fprintf(' Priors do not have right dimensions at restart. \n')
    fprintf(' Modify DefineInputsForInverseRun to ensure that dimensions are correct.\n')
    error('Ua:GetInputForInverseRestartRun:incorrectdimensions','incorrect dimensions')
    
end


% [UserVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint]=DefineModificationsToInverseRestartRunData(UserVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint) ; 



% %% TestIng
% BCsAdjoint.ubFixedNode=MUA.Boundary.Nodes ;   BCsAdjoint.ubFixedValue=BCsAdjoint.ubFixedNode*0;
% BCsAdjoint.vbFixedNode=MUA.Boundary.Nodes ;   BCsAdjoint.vbFixedValue=BCsAdjoint.vbFixedNode*0;
% figure ; PlotBoundaryConditions(CtrlVar,MUA,BCsAdjoint) ;
% %%







end