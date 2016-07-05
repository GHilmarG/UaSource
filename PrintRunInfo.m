function PrintRunInfo(CtrlVar)

fprintf(CtrlVar.fidlog,'\n \n');
fprintf(CtrlVar.fidlog,' **************************     Úa version: 5 July 2016  ********************** \n');
fprintf(CtrlVar.fidlog,'    Run starts at %s  \n ',datestr(now));
fprintf(CtrlVar.fidlog,'   Experiment: %s  \n ',CtrlVar.Experiment);

if CtrlVar.doInverseStep
    if CtrlVar.Restart
        fprintf(CtrlVar.fidlog,'   Inverse-modelling restart run.  \n ');
    else
        fprintf(CtrlVar.fidlog,'   Inverse-modelling run.  \n ');
    end
end


if  CtrlVar.doDiagnostic
    if CtrlVar.Restart
        fprintf(CtrlVar.fidlog,'   Time-independent restart run.  \n ');
    else
        fprintf(CtrlVar.fidlog,'   Time-independent model run.  \n ');
    end
end

if  CtrlVar.doPrognostic
    if CtrlVar.Restart
        fprintf(CtrlVar.fidlog,'   Time-dependent restart run.  \n ');
    else
        fprintf(CtrlVar.fidlog,'   Time-dependent model run.  \n ');
    end
end


[status,hostname]=system('hostname');
if status==0
    fprintf('   Running on: %s ',hostname)
end

fprintf(CtrlVar.fidlog,'***************************************************************************** \n \n');


end
