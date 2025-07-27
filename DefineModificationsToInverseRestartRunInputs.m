



function [UserVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo]=...
    DefineModificationsToInverseRestartRunInputs(UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo) 


%% This file can be used to modify input data for an inverse restart run.
%
% This is a bit unusual thing to do, and potentially dangerous. 
%
% But this is also a convenient way of changing some fields used in a restart inverse run.
%
%
% For example, if some of the measurement errors that were initially defined in 
%
%   DefineInputsForInverseRun.m
%
% the one could here modify some of the fields of 
%
%   Meas
%
%
% For example, to modify the data errors: 
% 
%   dhdtErr=1+zeros(MUA.Nnodes,1);
%   Meas.dhdtCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,dhdtErr.^2,MUA.Nnodes,MUA.Nnodes);
% 
%   Error=1+zeros(MUA.Nnodes,1);
%   Meas.usCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,Error.^2,MUA.Nnodes,MUA.Nnodes);
%   Meas.vsCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,Error.^2,MUA.Nnodes,MUA.Nnodes);
% 
%%



end
