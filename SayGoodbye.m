function SayGoodbye(CtrlVar,RunInfo)
    
    if CtrlVar.WriteRunInfoFile
        fclose(RunInfo.File.fid);
    end
    fprintf(CtrlVar.fidlog,' Run finishes at %s \n ================           Allt gott �� endirinn er allra bestur!    ======================\n \n',datetime);
    
end