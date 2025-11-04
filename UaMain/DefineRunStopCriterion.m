function [UserVar,isStop]=DefineRunStopCriterion(UserVar,RunInfo,CtrlVar,MUA,BCs,F) 

%% DefineRunStopCriterion.m
%
% If the field 
%
%   CtrlVar.UseUserDefinedRunStopCriterion 
%
% is set to true in 
%
%   Ua2D_InitialUserInput.m
% 
% then the m-file 
%
%   DefineRunStopCriterion.m
%
% can be used to specify a run-stop criterion. Make sure this m-file is in your
% local run directory. 
%
% The criterion is used within the run-step loop, ie in both time dependent and
% time independent runs.
%
% Note that this criterion is not used in an inverse run and has therefore no
% effect on if and when an inverse run stops.
%
%%

isStop=false;

%% Example of a run-stop criterion based on norm of dh/dt.
%
%   RSC=norm(F.dhdt)/sqrt(numel(F.dhdt));
%   fprintf(' norm(F.dhdt)/sqrt(numel(F.dhdt)) = %f \n ',RSC)
% 
%   if RSC==0 || isempty(RSC)  % If this is exactly zero, then dhdt has presumably not been calculated 
%     isStop=false; 
%   else
%     if RSC < 10 
%         isStop=true;
%     end
%   end



end