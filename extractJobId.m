function jobID = extractJobId(cmdOut)
% Extracts the job ID from the qsub command output for SGE

% Copyright 2010 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2010/03/22 03:43:18 $

% The output of qsub will be:
% Your job 496 ("<job name>") has been submitted

% Now parse the output of bsub to extract the job number
jobNumberStr = regexp(cmdOut, 'job [0-9]*', 'once', 'match');
jobID = sscanf(jobNumberStr, 'job %d');
dctSchedulerMessage(0, '%s: Job ID %d was extracted from qsub output %s.', mfilename, jobID, cmdOut);