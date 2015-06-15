function state = getJobStateFcn(scheduler, job, state)
%GETJOBSTATEFCN Gets the state of a job from SGE
%
% Set your scheduler's GetJobStateFcn to this function using the following
% command:
%     set(sched, 'GetJobStateFcn', @getJobStateFcn);

% Copyright 2010 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2010/03/31 18:14:17 $

% Store the current filename for the dctSchedulerMessages
currFilename = mfilename;
if ~scheduler.HasSharedFilesystem
    error('distcompexamples:SGE:SubmitFcnError', ...
        'The submit function %s is for use with shared filesystems.', currFilename)
end


% Shortcut if the job state is already finished or failed
jobInTerminalState = strcmp(state, 'finished') || strcmp(state, 'failed');
if jobInTerminalState
    return;
end
 % Get the information about the actual scheduler used
data = scheduler.getJobSchedulerData(job);
if isempty(data)
    % This indicates that the job has not been submitted, so just return
    dctSchedulerMessage(1, '%s: Job scheduler data was empty for job with ID %d.', currFilename, job.ID);
    return
end
try
    jobIDs = data.SchedulerJobIDs;
catch err
    ex = MException('distcompexamples:SGE:FailedToRetrieveJobID', ...
        'Failed to retrieve scheduler''s job IDs from the job scheduler data.');
    ex = ex.addCause(err);
    throw(ex);
end
  
% Get the full xml from qstat so that we can look for 
% <job_list state="pending">
% <job_list state="running">
commandToRun = sprintf('qstat -xml');
dctSchedulerMessage(4, '%s: Querying scheduler for job state using command:\n\t%s', currFilename, commandToRun);

try
    % We will ignore the status returned from the state command because
    % a non-zero status is returned if the job no longer exists
    % Make the shelled out call to run the command.
    [~, cmdOut] = system(commandToRun);
catch err
    ex = MException('distcompexamples:SGE:FailedToGetJobState', ...
        'Failed to get job state from scheduler.');
    ex.addCause(err);
    throw(ex);
end

schedulerState = iExtractJobState(cmdOut, jobIDs);
dctSchedulerMessage(6, '%s: State %s was extracted from scheduler output:\n\t%s', currFilename, schedulerState, cmdOut);

% If we could determine the scheduler's state, we'll use that, otherwise
% stick with MATLAB's job state.
if ~strcmp(schedulerState, 'unknown')
    state = schedulerState;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = iExtractJobState(qstatXMLOut, requestedJobIDs)
% Function to extract the job state for the requested jobs from the 
% output of qstat -xml

% write the xml to a temp file first so that we can use xmlread
filename = tempname;
fid = fopen(filename, 'w');
% If we couldn't write to the file, then just return state unknown
if fid < 0
    state = 'unknown';
    return;
end
fprintf(fid, qstatXMLOut);
fclose(fid);

% use xmlread to read in the file
xmlFileDOM = xmlread(filename);
% now delete the temporary file
delete(filename);

% We expect the XML to be in the following format:
%  <queue_info>
%    <job_list state="running">
%      <JB_job_number>535</JB_job_number>
%      <JAT_prio>0.55500</JAT_prio>
%      <JB_name>Job1.1</JB_name>
%      <JB_owner>elwinc</JB_owner>
%      <state>r</state>
%      <JAT_start_time>2010-01-28T05:56:27</JAT_start_time>
%      <queue_name>all.q@dct13glnxa64.mathworks.com</queue_name>
%      <slots>1</slots>
%    </job_list>
%
% Find the correct JB_job_number node and get the parent to find out 
% the job's state
% the job's state
allJobNodes = xmlFileDOM.getElementsByTagName('JB_job_number');

numJobsFound = allJobNodes.getLength();
if numJobsFound == 0
    % No jobs in the qstat output, so the one we are interested in
    % must be finished
    state = 'finished';
    return;
end

stateParseError = MException('distcompexamples:SGE:StateParseError', ...
    'Failed to parse XML output from qstat.  Could not find "state" attribute for job_list node.');

for ii = 0:allJobNodes.getLength() - 1
    jobNode = allJobNodes.item(ii);
    jobNumber = str2double(jobNode.getFirstChild.getData);
    
    isEqualFcn = @(x) isequal(jobNumber, x);
    if any(cellfun(isEqualFcn, requestedJobIDs))
        % Only get the parent node if the current node is a
        % job that we are interested in.
        jobParentNode = jobNode.getParentNode;
        % We were expecting the state attribute
        if ~jobParentNode.hasAttributes
            throw(stateParseError);
        end
        
        stateAttribute = jobParentNode.getAttributes.getNamedItem('state');
        if isempty(stateAttribute)
            throw(stateParseError);
        end
        
        % We've got the state for the job that we are interested in.
        jobState = char(stateAttribute.getValue);
        if strcmpi(jobState, 'running');
            % If any of the requested jobs are running, then the whole
            % job is still running, so we can stop searching and just
            %return with state running
            state = 'running';
            return;
        end
    end
end

state = 'unknown';
