function isSelected=uvhNewtonStateSelected(Opt,CtrlVar,iteration)

%% Is the current uvh Newton-Raphson iteration selected for saving or testing?
%
%   isSelected=uvhNewtonStateSelected(Opt,CtrlVar,iteration)
%
% (15 Oct 2026) Used for CtrlVar.uvhSaveNewtonSystem and CtrlVar.uvhTestJacobian (see Ua2D_DefaultParameters.m). An iteration is selected
% if the option is on, the Newton-Raphson iteration number is in Opt.Iterations, and the run step is in Opt.RunSteps. An empty
% Opt.Iterations or Opt.RunSteps means all.
%
%%

isSelected=false ;

if isempty(Opt) || ~isstruct(Opt)
    return
end

if isfield(Opt,"Save")
    isOn=Opt.Save ;
elseif isfield(Opt,"Test")
    isOn=Opt.Test ;
else
    isOn=false ;
end

if ~isOn
    return
end

if isfield(Opt,"Iterations") && ~isempty(Opt.Iterations) && ~ismember(iteration,Opt.Iterations)
    return
end

if isfield(Opt,"RunSteps") && ~isempty(Opt.RunSteps) && ~ismember(CtrlVar.CurrentRunStepNumber,Opt.RunSteps)
    return
end

isSelected=true ;

end