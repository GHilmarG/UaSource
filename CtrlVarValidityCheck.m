function CtrlVar=CtrlVarValidityCheck(CtrlVar)

%  Performs some basic validity checks on CtrlVar

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
end

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
    
    if strcmpi(CtrlVar.Inverse.DataMisfit.GradientCalculation,'fixpoint')
        
        % if fixpoint, then only c inversion is possible
        CtrlVar.Inverse.Regularize.Field=replace(lower(CtrlVar.Inverse.Regularize.Field),'logaglen','');
        CtrlVar.Inverse.Regularize.Field=replace(lower(CtrlVar.Inverse.Regularize.Field),'aglen','');
        CtrlVar.Inverse.InvertFor=replace(lower(CtrlVar.Inverse.InvertFor),'logaglen','');
        CtrlVar.Inverse.InvertFor=replace(lower(CtrlVar.Inverse.InvertFor),'aglen','');
        
    end
    
    if ~contains(lower(CtrlVar.Inverse.InvertFor),'aglen')
        
        CtrlVar.Inverse.Regularize.Field=replace(lower(CtrlVar.Inverse.Regularize.Field),'logaglen','');
        CtrlVar.Inverse.Regularize.Field=replace(lower(CtrlVar.Inverse.Regularize.Field),'aglen','');
        
    end
    
    if ~contains(lower(CtrlVar.Inverse.InvertFor),'c')
        
        CtrlVar.Inverse.Regularize.Field=replace(lower(CtrlVar.Inverse.Regularize.Field),'logc','');
        CtrlVar.Inverse.Regularize.Field=replace(lower(CtrlVar.Inverse.Regularize.Field),'c','');
        
    end
    
 
    
    if ~contains(lower(CtrlVar.Inverse.InvertFor),["aglen","c"])
        fprintf('the string CtrlVar.Inverse.InvertFor must contain ``AGlen`` and/or ``C`` \n')
        error('Invalid inputs.')
    end
    
    if ~contains(lower(CtrlVar.Inverse.DataMisfit.GradientCalculation),["fixpoint","adjoint"])
        fprintf('the string CtrlVar.Inverse.DataMisfit.GradientCalculation must contain either ''fixpoint'' or ``adjoint` \n')
        error('Invalid inputs.')
    end
    
    
end


if isfield(CtrlVar,'AdaptMeshInterval')
    
    fprintf(' Note: CtrlVar.AdaptMeshInterval no longer used. Use CtrlVar.AdaptMeshRunStepInterval instead.\n')
    fprintf(' setting CtrlVar.AdaptMeshRunStepInterval=CtrlVar.AdaptMeshInterval \n')
    CtrlVar.AdaptMeshRunStepInterval=CtrlVar.AdaptMeshInterval;
    
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





end

