
function WriteForwardRunRestartFile(UserVar,CtrlVar,MUA,BCs,F,GF,l,RunInfo)

RestartFile=CtrlVar.NameOfRestartFiletoWrite;
fprintf(CtrlVar.fidlog,' \n ################## %s %s ################### \n Writing restart file %s  at t=%-g \n %s \n ',CtrlVar.Experiment,datestr(now),RestartFile,CtrlVar.time);



CtrlVarInRestartFile=CtrlVar;
UserVarInRestartFile=UserVar;

if isfield(RunInfo.Forward,'dtRestart')
    CtrlVarInRestartFile.dt=RunInfo.Forward.dtRestart;
end


try
    save(RestartFile,'CtrlVarInRestartFile','UserVarInRestartFile','MUA','BCs','F','GF','l','RunInfo','-v7.3');
    fprintf(CtrlVar.fidlog,' Writing restart file was successful. \n');

catch exception

    warning("WriteForwardRunRestartFile:RestartFileNotSaved","Could not save restart file %s \n ",RestartFile)
    fprintf(CtrlVar.fidlog,' Could not save restart file %s \n ',RestartFile);
    fprintf(CtrlVar.fidlog,'%s \n',exception.message);

    try
        RestartFile="ForwardRestartFile.mat";
        save(RestartFile,'CtrlVarInRestartFile','UserVarInRestartFile','MUA','BCs','F','GF','l','RunInfo','-v7.3');
        fprintf(CtrlVar.fidlog,' Writing local restart file %s was successful. \n',RestartFile);

    catch exception

        % OK, try to write a default restart file in local directory
        warning("WriteForwardRunRestartFile:RestartFileNotSaved","Could not save restart file %s \n ",RestartFile)
        fprintf(CtrlVar.fidlog,' Could not save restart file %s \n ',RestartFile);
        fprintf(CtrlVar.fidlog,'%s \n',exception.message);

    end


end


end
