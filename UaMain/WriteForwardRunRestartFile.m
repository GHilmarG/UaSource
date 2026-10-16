




function WriteForwardRunRestartFile(UserVar,CtrlVar,MUA,BCs,F,GF,l,RunInfo,Fm1)

% (9 Oct 2026) Fm1, the rates of the previous time step (see UpdateFtimeDerivatives.m), is also saved. It is needed for the
% second-order explicit estimate (ExplicitEstimationUsingBackwardDifferences.m) in the first time step after a restart.
if nargin<9 || isempty(Fm1)
    Fm1=UaFields;
end

RestartFile=CtrlVar.NameOfRestartFiletoWrite;

if isempty(RestartFile) || RestartFile==""
    return
end


fprintf("\n ################## Experiment: %s ################### \n Writing restart file: %s  \n for simulation time t=%-g at wall-clock time %s\n \n ",CtrlVar.Experiment,CtrlVar.NameOfRestartFiletoWrite,CtrlVar.time,datetime);


MUA.workers=[]; % saving composites is not supported, MATLAB 2024

CtrlVarInRestartFile=CtrlVar;
UserVarInRestartFile=UserVar;

if isfield(RunInfo.Forward,'dtRestart')
    CtrlVarInRestartFile.dt=RunInfo.Forward.dtRestart;
end


try
    
    save(RestartFile,'CtrlVarInRestartFile','UserVarInRestartFile','MUA','BCs','F','Fm1','GF','l','RunInfo','-v7.3');
    fprintf(CtrlVar.fidlog,' Writing restart file was successful. \n');

catch exception

    warning("WriteForwardRunRestartFile:RestartFileNotSaved","Could not save restart file %s \n ",RestartFile)
    fprintf(CtrlVar.fidlog,' Could not save restart file %s \n ',RestartFile);
    fprintf(CtrlVar.fidlog,'%s \n',exception.message);

    try
        RestartFile="ForwardRestartFile.mat";
        save(RestartFile,'CtrlVarInRestartFile','UserVarInRestartFile','MUA','BCs','F','Fm1','GF','l','RunInfo','-v7.3');
        fprintf(CtrlVar.fidlog,' Writing local restart file %s was successful. \n',RestartFile);

    catch exception

        % OK, try to write a default restart file in local directory
        warning("WriteForwardRunRestartFile:RestartFileNotSaved","Could not save restart file %s \n ",RestartFile)
        fprintf(CtrlVar.fidlog,' Could not save restart file %s \n ',RestartFile);
        fprintf(CtrlVar.fidlog,'%s \n',exception.message);

    end


end


end
