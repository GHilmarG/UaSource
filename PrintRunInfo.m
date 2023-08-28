
function PrintRunInfo(CtrlVar)


fprintf(CtrlVar.fidlog,'\n \n');
fprintf(CtrlVar.fidlog,' **************************    Úa version: 28 August 2023 (beta)      ********************** \n');
fprintf(CtrlVar.fidlog,'    Run starts at %s  \n ',datetime);
fprintf(CtrlVar.fidlog,'   Experiment: %s  \n ',CtrlVar.Experiment);

if CtrlVar.doInverseStep
    if CtrlVar.Restart
        fprintf(CtrlVar.fidlog,'   Inverse-modelling restart run.  \n ');
    else
        fprintf(CtrlVar.fidlog,'   Inverse-modelling run.  \n ');
    end
    
    fprintf('\tCtrlVar.Inverse.MinimisationMethod=''%s''\n' ,CtrlVar.Inverse.MinimisationMethod)
    fprintf('\tCtrlVar.Inverse.InvertFor=''%s''\n',CtrlVar.Inverse.InvertFor)
    fprintf('\tCtrlVar.Inverse.Measurements=''%s''\n',CtrlVar.Inverse.Measurements)
    fprintf('\tCtrlVar.Inverse.Regularize.Field=''%s''\n',CtrlVar.Inverse.Regularize.Field)
    fprintf('\tCtrlVar.Inverse.DataMisfit.GradientCalculation=''%s''\n',CtrlVar.Inverse.DataMisfit.GradientCalculation)
    
    
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
if  status==0
    fprintf('   Running on: %s ',hostname)
end

fprintf(CtrlVar.fidlog,'***************************************************************************** \n \n');


end
