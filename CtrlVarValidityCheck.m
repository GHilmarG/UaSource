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

% if  a restart and an inverse run
% then interpret this as a restarted inversion and
% set AdjointRestart to true and Restart to false
if CtrlVar.Restart && CtrlVar.InverseRun
    CtrlVar.AdjointRestart=1;
    CtrlVar.Restart=0;
end


if CtrlVar.AdjointRestart && ~CtrlVar.InverseRun
    fprintf('CtrlVar.AdjointRestart set to true, but CtrlVar.InverseRun not.\n')
    fprintf('Set either AdjointRestart to false or InverseRun to true.\n')
    error('Inconsistent definition of CtrlVar')
end

if ~ismember(CtrlVar.TriNodes,[3 6 10])
    fprintf(CtrlVar.fidlog,'CtrlVar.Trinodes=%-i but must be either 3, 6 or 10 \n',CtrlVar.TriNodes);
    error('CtrlVarValidityCheck:Trinodes','CtrlVar not valid')
end

if ~iscell(CtrlVar.RefineCriteria) ; CtrlVar.RefineCriteria=cellstr(CtrlVar.RefineCriteria) ; end

[n1,n2]=size(CtrlVar.RefineCriteria);
if ~isfield(CtrlVar,'RefineCriteriaWeights')
    CtrlVar.RefineCriteriaWeights=zeros(1,n2)+1;
end



if ~isfield(CtrlVar,'RefineCriteriaFlotationLimit')
    CtrlVar.RefineCriteriaFlotationLimit=zeros(1,n2)+NaN  ;
end


[m1,m2]=size(CtrlVar.RefineCriteriaFlotationLimit);

if n2~=m2 ;
    fprintf(CtrlVar.fidlog,'Number of elements in CtrlVar.RefineCriteriaFloationLimit (%i) not equal to number of RefineCriteria (%i) \n',m2,n2);
    if n2==1
        fprintf(CtrlVar.fidlog,'Assuming no such limits for any of the mesh refinement criteria. \n');
        CtrlVar.RefineCriteriaFlotationLimit=zeros(n2,1)+NaN;
    else
        error('CtrlVarValidityCheck:RefineCriteriaFloationLimit','CtrlVar not valid')
    end
end

[m1,m2]=size(CtrlVar.RefineCriteriaWeights);
if n2~=m2 ;
    fprintf(CtrlVar.fidlog,'Number of elements in CtrlVar.RefineCriteriaWeights (%i) not equal to number of RefineCriteria (%i) \n',m2,n2);
    if n2==1
        CtrlVar.RefineCriteriaWeights=zeros(n2,1)+1;
        fprintf(' Assuming all criteria have equal weights. All weights set to 1 \n')
    else
        error('CtrlVarValidityCheck:RefineCriteriaWeights','CtrlVar not valid')
    end
end

%% if user requires an inverse step, set prognostic and diagnostic flags to zero
if CtrlVar.InverseRun==1    % inverse step takes precedence over prognostic and diagnostic, if conflict
    CtrlVar.doDiagnostic=0  ;
    CtrlVar.doPrognostic=0 ;
    CtrlVar.AdaptMesh=0;
    CtrlVar.Restart=0;
else
    CtrlVar.AdjointRestart=0 ;
    if CtrlVar.doDiagnostic  ; % diagnostic takes precedence over prognostic, if conflict
        CtrlVar.doPrognostic=0 ;
        CtrlVar.nTimeSteps=1;
    end
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

end

