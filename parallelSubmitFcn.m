function parallelSubmitFcn(scheduler, job, props)
%PARALLELSUBMITFCN Submit a parallel MATLAB job to a SGE scheduler
%
% Set your scheduler's ParallelSubmitFcn to this function using the following
% command:
%     set(sched, 'ParallelSubmitFcn', @parallelSubmitFcn);
%
% See also parallel.cluster.generic.parallelDecodeFcn.
%

% Copyright 2010 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2010/03/31 18:14:18 $

decodeFunction = 'parallel.cluster.generic.parallelDecodeFcn';

% Store the current filename for the dctSchedulerMessages
currFilename = mfilename;
if ~scheduler.HasSharedFilesystem
    error('distcompexamples:SGE:SubmitFcnError', ...
        'The submit function %s is for use with shared filesystems.', currFilename)
end

% Ensure that the cluster size is consistent with the job's 
% minimum number of workers
minProcessors = job.MinimumNumberOfWorkers;
if minProcessors > scheduler.ClusterSize
    error('distcompexamples:SGE:ResourceLimit', ...
        ['You requested a minimum of %d workers, but the scheduler''s ClusterSize property ' ...
        'is currently set to allow a maximum of %d workers.  ' ...
        'To run a parallel job with more tasks than this, increase the value of the ClusterSize ' ...
        'property for the scheduler.'], ...
        minProcessors, scheduler.ClusterSize);
end

if ~strcmpi(scheduler.ClusterOsType, 'unix')
    error('distcompexamples:SGE:SubmitFcnError', ...
        'The submit function %s only supports clusters with unix OS.', currFilename)
end

% The job specific environment variables
% Remove leading and trailing whitespace from the MATLAB arguments
matlabArguments = strtrim(props.MatlabArguments);
variables = {'MDCE_DECODE_FUNCTION', decodeFunction; ...
    'MDCE_STORAGE_CONSTRUCTOR', props.StorageConstructor; ...
    'MDCE_JOB_LOCATION', props.JobLocation; ...
    'MDCE_MATLAB_EXE', props.MatlabExecutable; ... 
    'MDCE_MATLAB_ARGS', matlabArguments; ...
    'MDCE_DEBUG', 'true'; ...
    'MDCE_STORAGE_LOCATION', props.StorageLocation; ...
    'MDCE_CMR', scheduler.ClusterMatlabRoot; ...
    'MDCE_TOTAL_TASKS', num2str(props.NumberOfTasks)};
% Set the required environment variables
for ii = 1:size(variables, 1)
    setenv(variables{ii,1}, variables{ii,2});
end
% Which variables do we need to forward for the job?  Take all
% those that we have set.
variablesToForward = variables(:,1);

% Deduce the correct quote to use based on the OS of the current machine
if ispc
    quote = '"';
else 
    quote = '''';
end

% The local job directory
localJobDirectory = fullfile(scheduler.DataLocation, props.JobLocation);


% The script name is parallelJobWrapper.sh
scriptName = 'parallelJobWrapper.sh';
% The wrapper script is in the same directory as this file
dirpart = fileparts(mfilename('fullpath'));
quotedScriptName = sprintf('%s%s%s', quote, fullfile(dirpart, scriptName), quote);

% Choose a file for the output. Please note that currently, DataLocation refers
% to a directory on disk, but this may change in the future.
logFile = fullfile(localJobDirectory, sprintf('Job%d.mpiexec.log', job.ID));
quotedLogFile = sprintf('%s%s%s', quote, logFile, quote);

jobName = sprintf('Job%d', job.ID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CUSTOMIZATION MAY BE REQUIRED %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose a number of slots per node to use 
% You may wish to customize this section to match your cluster, 
% for example if you wish to limit the number of slots that 
% can be used for a single job.
procsPerNode = 1;
numberOfNodes = ceil(props.NumberOfTasks/procsPerNode);
% You may also with to supply additional submission arguments to 
% the qsub command here.
%additionalSubmitArgs  = sprintf('-q compute.q -pe matlab_2011b %d', numberOfNodes);
%
% BAS Customisation.
%
matlabQueue = lower(getenv('BAS_MATLAB_QUEUE'))
if not (strcmp(matlabQueue,'compute.q')||strcmp(matlabQueue,'cuda.q'))
  matlabQueue = 'compute.q';
end
additionalSubmitArgs  = sprintf('-q %s -pe matlab_2011b %d -M %s@bas.ac.uk -m bea', matlabQueue,numberOfNodes,getenv('USER'));
%
% End BAS custinisation.
%

dctSchedulerMessage(4, '%s: Requesting %d nodes with %d processors per node', currFilename, ...
    numberOfNodes, procsPerNode);
dctSchedulerMessage(5, '%s: Generating command for task %i', currFilename, ii);
commandToRun = getSubmitString(jobName, quotedLogFile, quotedScriptName, ...
    variablesToForward, additionalSubmitArgs);   


% Now ask the cluster to run the submission command
dctSchedulerMessage(4, '%s: Submitting job using command:\n\t%s', currFilename, commandToRun);
try
    % Make the shelled out call to run the command.
    [cmdFailed, cmdOut] = system(commandToRun);
catch err
    cmdFailed = true;
    cmdOut = err.message;
end
if cmdFailed
    error('distcompexamples:SGE:SubmissionFailed', ...
        'Submit failed with the following message:\n%s', cmdOut);
end

dctSchedulerMessage(1, '%s: Job output will be written to: %s\nSubmission output: %s\n', currFilename, logFile, cmdOut);

jobIDs = extractJobId(cmdOut);
% jobIDs must be a cell array
if isempty(jobIDs)
    warning('distcompexamples:SGE:FailedToParseSubmissionOutput', ...
        'Failed to parse the job identifier from the submission output: "%s"', ...
        cmdOut);
end
if ~iscell(jobIDs)
    jobIDs = {jobIDs};
end

% set the job ID on the job scheduler data
scheduler.setJobSchedulerData(job, struct('SchedulerJobIDs', {jobIDs}));


