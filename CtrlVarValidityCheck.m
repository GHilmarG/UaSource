function CtrlVar=CtrlVarValidityCheck(CtrlVar)

%  Performs some basic validity checks on CtrlVar


fn = fieldnames(CtrlVar.MustBe);
for I=1:numel(fn)
    CheckUaCtrlVarFields(CtrlVar,fn{I})
end

%% is user still using doDiagnosti/doPrognostic, set TimeDependentRun accordingly
if isfield(CtrlVar,'doDiagnostic')
    if CtrlVar.doDiagnostic
        CtrlVar.TimeDependentRun=0 ;
    else
        CtrlVar.TimeDependentRun=1 ;
    end
end


if isfield(CtrlVar,'doPrognostic')
    if CtrlVar.doPrognostic
        CtrlVar.TimeDependentRun=1 ;
    else
        CtrlVar.TimeDependentRun=0 ;
    end
end

if isfield(CtrlVar,'doInverseStep')
    if CtrlVar.doInverseStep
        CtrlVar.InverseRun=1;
    else
        CtrlVar.InverseRun=0;
    end
end




%%

if CtrlVar.TimeDependentRun
    CtrlVar.doDiagnostic=0  ;
    CtrlVar.doPrognostic=1 ;
else
    CtrlVar.doDiagnostic=1  ;
    CtrlVar.doPrognostic=0 ;
end

if CtrlVar.InverseRun
    CtrlVar.doInverseStep=1 ;
else
    CtrlVar.doInverseStep=0 ;
end

if ~ismember(CtrlVar.TriNodes,[3 6 10])
    fprintf(CtrlVar.fidlog,'CtrlVar.Trinodes=%-i but must be either 3, 6 or 10 \n',CtrlVar.TriNodes);
    error('Ua:CtrlVarValidityCheck:Trinodes','CtrlVar not valid')
end


%% if user requires an inverse step, set prognostic and diagnostic flags to zero
if CtrlVar.InverseRun    % inverse step takes precedence over prognostic and diagnostic, if conflict
    CtrlVar.doDiagnostic=0  ;
    CtrlVar.doPrognostic=0 ;
    CtrlVar.AdaptMesh=0;
    CtrlVar.MUA.MassMatrix=true;
    CtrlVar.MUA.StiffnessMatrix=true;
end


if CtrlVar.InverseRun  && contains(CtrlVar.SlidingLaw,["Coulomb","-C-"])
    fprintf("Inversion using %s sliding law is not possible!\n",CtrlVar.SlidingLaw)
    error('Ua:CtrlVarValidityCheck:InverseAdapt','CtrlVar not valid')
end

%%

if CtrlVar.UaRunType==""  %  "-uvh-" , "-uv-h-" ,  "-uv-" , "-h-" ;

    if  CtrlVar.TimeDependentRun

        if CtrlVar.Implicituvh
            CtrlVar.UaRunType="-uvh-";
        else
            CtrlVar.UaRunType="-uv-h-";
        end

    else

        CtrlVar.UaRunType="-uv-";

    end

else

    if CtrlVar.UaRunType=="-uv-"
        CtrlVar.TimeDependentRun=0;
    else
        CtrlVar.TimeDependentRun=1;

        if CtrlVar.UaRunType=="-uvh-"
            CtrlVar.Implicituvh=1;
        else
            CtrlVar.Implicituvh=0;
        end
    end

end

%%




if CtrlVar.doplots==0 ; CtrlVar.PlotMesh=0; end

if CtrlVar.StandartOutToLogfile
    CtrlVar.fidlog=fopen(CtrlVar.Logfile,'w');
else
    CtrlVar.fidlog=1;
end

% make sure that Cmin is zero in the SSHEET and the Hybrid flow approximation
if strcmpi(CtrlVar.FlowApproximation,'ssheet') || strcmpi(CtrlVar.FlowApproximation,'hybrid')
    if CtrlVar.Cmin>0
        fprintf('In the SSHEET and the Hybrid flow approximations Cmin must be zero. Cmin set to zero.\n')
        CtrlVar.Cmin=0;
    end
end

if CtrlVar.TimeDependentRun &&  ~CtrlVar.Restart
    CtrlVar.InitialDiagnosticStep=1;  % Always do an initial diagnostic step when starting a transient run. But if a restart run, then do not change this value from that set by the user.
    % An inital diagnostic step is therefore done if:
    % 1) so asked by the user, ie if the user sets CtrlVar.InitialDiagnosticStep=1, and
    % 2) at the beginning of an implicut uvh transient run.
    % Unless asked by the user, no initial diagnostic step is done at the beginning of an implicut uvh transient restart run.
end

% AdaptMesh not allowed in combination with an inverse run

if CtrlVar.AdaptMesh && CtrlVar.InverseRun
    fprintf('UaError: Both CtrlVar.AdaptMesh and CtrlVar.InverseRun are set to true!\n')
    fprintf('As adaptive meshing can only be done in a combination with a forward run this combination is not allowed.\')
    error('Ua:CtrlVarValidityCheck:InverseAdapt','CtrlVar not valid')
    
end


if isfield(CtrlVar,'nTimeSteps')
    warning('Ua:CtrlVarValidityCheck:nTimeSteps','The field nTimeSteps in the strucure array CtlrVar no longer used. Use the TotalNumberOfForwardRunSteps field instead')
    CtrlVar.TotalNumberOfForwardRunSteps=CtrlVar.nTimeSteps;
end

if isfield(CtrlVar,'GmeshFile')
    error('Ua:CtrlVarValidyCheck','The field CtrlVar.GmeshFile no longer used. Replace with CtrlVar.GmshFile.')
end


if isfield(CtrlVar,'GmeshMeshingAlgorithm')
    error('Ua:CtrlVarValidyCheck','The field CtrlVar.GmeshMeshingAlgorithm no longer used. Replace with CtrlVar.GmshMeshingAlgorithm.')
end

%% adapt

if isfield(CtrlVar,'AdaptMeshIterations')
    
    fprintf('Note: CtrlVar.AdaptMeshIterations is no longer used.\n')
    fprintf('Use   CtrlVar.AdaptMeshMaxIterations instead.\n')
    error('Ua:CtrlVarValidyCheck','The field CtrlVar.AdaptMeshIterations is no longer used. Replace with CtrlVar.AdaptMeshMaxIterations')
    
end





%% inverse

if CtrlVar.InverseRun
    
    % First make sure that CtrlVar.Inverse.InvertFor and CtrlVar.Inverse.Regularize.Field
    % only contain some combinations of "-C-","-logC-","-AGlen-","-logAGlen-" and
    % "-B-"
    %
    % for example : "-A-C-"  -> "-AGlen-C-"
    
    [CtrlVar.Inverse.InvertFor,status]=SearchAndReplaceInverseFieldsInCtrlVar(CtrlVar.Inverse.InvertFor);
    
    if ~status
        fprintf(" CtrlVar.Inverse.InvertFor does not appear to have a valid value.\n")
        fprintf(" CtrlVar.Inverse.InvertFor=%s \n",CtrlVar.Inverse.InvertFor)
        error("CtrlVarValidityCheck:CtrlVar.Inverse.InvertForInvalid")
    end
    
    
    
    [CtrlVar.Inverse.Regularize.Field,status]=SearchAndReplaceInverseFieldsInCtrlVar(CtrlVar.Inverse.Regularize.Field);
    
    if ~status
        fprintf(" CtrlVar.Inverse.Regularize.Field does not appear to have a valid value.\n")
        fprintf(" CtrlVar.Inverse.Regularize.Field=%s \n",CtrlVar.Inverse.Regularize.Field);
        error("CtrlVarValidityCheck:CtrlVar.Inverse.Regularize.Field")
    end
    
    if strcmpi(CtrlVar.Inverse.DataMisfit.GradientCalculation,"fixpoint")
        
        % if fixpoint, then only c inversion is possible
        CtrlVar.Inverse.Regularize.Field=replace(CtrlVar.Inverse.Regularize.Field,"logAGlen","");
        CtrlVar.Inverse.Regularize.Field=replace(CtrlVar.Inverse.Regularize.Field,"Aglen","");
        CtrlVar.Inverse.InvertFor=replace(CtrlVar.Inverse.InvertFor,"logAGlen","");
        CtrlVar.Inverse.InvertFor=replace(CtrlVar.Inverse.InvertFor,"AGlen","");
        
        
        
    end
    
    % Don't regularize A if not inverting for A, so
    if ~contains(CtrlVar.Inverse.InvertFor,"AGlen")
        
        CtrlVar.Inverse.Regularize.Field=replace(CtrlVar.Inverse.Regularize.Field,"-logAGlen-","");
        CtrlVar.Inverse.Regularize.Field=replace(CtrlVar.Inverse.Regularize.Field,"-AGlen-","");
        
        
    end
    
    % Don't regularize C if not inverting for C
    if ~contains(CtrlVar.Inverse.InvertFor,"C")
        
        CtrlVar.Inverse.Regularize.Field=replace(CtrlVar.Inverse.Regularize.Field,"-logC-","");
        CtrlVar.Inverse.Regularize.Field=replace(CtrlVar.Inverse.Regularize.Field,"-C-","");


    end

    CtrlVar.Inverse.Regularize.Field=replace(CtrlVar.Inverse.Regularize.Field,"--","-");
    CtrlVar.Inverse.InvertFor=replace(CtrlVar.Inverse.InvertFor,"--","-");

    if ~contains(lower(CtrlVar.Inverse.DataMisfit.GradientCalculation),["fixpoint","adjoint"])
        fprintf('the string CtrlVar.Inverse.DataMisfit.GradientCalculation must contain either ''fixpoint'' or ``adjoint` \n')
        error('Invalid inputs.')
    end

    % create a string with letters indicating which fields are being inverted for
    % e.g "ABC" if inverting for AGlen, B and C.
    CtrlVar.Inverse.InvertForField=string(sort(char(replace(replace(replace(string(CtrlVar.Inverse.InvertFor),"log","") ,"-",""),"AGlen","A")))) ;


    if contains(CtrlVar.Inverse.MinimisationMethod,"MatlabOptimization")

        if CtrlVar.Inverse.MinimisationMethod=="MatlabOptimization"
            fprintf("Inversion is HessianBased, ie provides a Hessian approximation. \n")
            CtrlVar.Inverse.MinimisationMethod="MatlabOptimization-HessianBased";
        end

    end



    if contains(CtrlVar.Inverse.MinimisationMethod,'Hessian')

        CtrlVar.Inverse.AdjointGradientPreMultiplier='I';

    end

end



if isfield(CtrlVar,'AdaptMeshInterval')

    fprintf(' Note: CtrlVar.AdaptMeshInterval no longer used. Use CtrlVar.AdaptMeshRunStepInterval instead.\n')
    error('Ua:CtrlVarValidityCheck','CtrlVar not valid')
    
end


if isfield(CtrlVar,'RefineCriteria')
    
    fprintf(' Note: CtrlVar.RefineCriteria no longer used. Use CtrlVar.ExplicitMeshRefinementCriteria instead.\n')
    error('Ua:CtrlVarValidityCheck','CtrlVar not valid')
    
end


if isfield(CtrlVar,'RefineCriteriaWeights')
    
    fprintf(' Note: CtrlVar.RefineCriteriaWeights no longer used. Use CtrlVar.ExplicitMeshRefinementCriteria instead.\n')
    error('Ua:CtrlVarValidityCheck','CtrlVar not valid')
    
end

if isfield(CtrlVar,'DefineOceanSurfaceAtEachTimeStep')
    
    fprintf(' Note: CtrlVar.DefineOceanSurfaceAtEachTimeStep no longer used.\n')
    fprintf('       Use CtrlVar.GeometricalVarsDefinedEachTransienRunStepByDefineGeometry instead.\n')
    fprintf('       For example: CtrlVar.GeometricalVarsDefinedEachTransienRunStepByDefineGeometry="S".\n')
    error('Ua:CtrlVarValidityCheck','CtrlVar not valid')
    
end


if isfield(CtrlVar,'InDiagnosticRunsDefineIceGeometryAtEveryRunStep')
    
    fprintf(' Note: CtrlVar.InDiagnosticRunsDefineIceGeometryAtEveryRunStep no longer used.\n')
    fprintf('       Use CtrlVar.GeometricalVarsDefinedEachDiagnosticRunStepByDefineGeometry instead.\n')
    fprintf('       For example: CtrlVar.GeometricalVarsDefinedEachDiagnosticRunStepByDefineGeometry="sbSB".\n')
    error('Ua:CtrlVarValidityCheck','CtrlVar not valid')
    
end

if isfield(CtrlVar,'ATStimeStepTarget')
    fprintf(' Note: CtrlVar.ATStimeStepTarget no longer used.\n')
    fprintf('       Use CtrlVar.ATSdtMax instead.\n')
    error('Ua:CtrlVarValidityCheck','CtrlVar not valid')
end

if isfield(CtrlVar,'UaOutputsDt')
    fprintf(' Note: CtrlVar.UaOutputsDt no longer used.\n')
    fprintf('       Use CtrlVar.DefineOutputsDt instead.\n')
    error('Ua:CtrlVarValidityCheck','CtrlVar not valid')
end

if isfield(CtrlVar.Inverse.TestAdjoint,'FiniteDifferenceType')
    fprintf(' Note: CtrlVar.Inverse.TestAdjoint.FiniteDifferenceType no longer used.\n')
    fprintf('       Use CtrlVar.TestAdjointFiniteDifferenceType instead.\n')
    error('Ua:CtrlVarValidityCheck','CtrlVar not valid')
end


if isfield(CtrlVar,"InfoLevelAdjoint")
    fprintf(' Note: CtrlVar.InfoLevelAdjoint no longer used.\n')
    fprintf('       Use CtrlVar.InfoLevelInverse instead.\n')
    error('Ua:CtrlVarValidityCheck','CtrlVar not valid')
end

if isfield(CtrlVar.Inverse,"MatlabOptimisationParameters")
    fprintf(' Note: CtrlVar.Inverse.MatlabOptimisationParameters no longer used.\n')
    fprintf("       Use \t CtrlVar.Inverse.MatlabOptimisationHessianParameters \n ")
    fprintf("        or \t CtrlVar.Inverse.MatlabOptimisationGradientParameters \n ")
    fprintf("        instead. \n")
    error('Ua:CtrlVarValidityCheck','CtrlVar not valid')
end



