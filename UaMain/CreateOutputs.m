



function [UserVar,RunInfo]=CreateOutputs(UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo)

narginchk(12,12)
nargoutchk(2,2)

warning('off','MATLAB:structOnObject')

% convert objects to structures. Otherwise it will not be possible to save/re-use these
% objects outside of �a.
%



RunInfo.WallTime.toc=datetime("now",Format="dd:hh:mm:ss.SSS");
RunInfo.WallTime.Total=duration(RunInfo.WallTime.Total+RunInfo.WallTime.toc-RunInfo.WallTime.tic,Format="dd:hh:mm:ss.SSS");
RunInfo.WallTime.tic=datetime("now",Format="dd:hh:mm:ss.SSS");

RunInfo.CPU.toc=duration(0,0,cputime,Format="dd:hh:mm:ss.SSS") ;
RunInfo.CPU.Total= RunInfo.CPU.Total+RunInfo.CPU.toc-RunInfo.CPU.tic; 
RunInfo.CPU.tic=duration(0,0,cputime,Format="dd:hh:mm:ss.SSS") ;

l=struct(l);
BCs=struct(BCs) ;
F=struct(F);
InvStartValues=struct(InvStartValues);
InvFinalValues=struct(InvFinalValues);
Priors=struct(Priors);
Meas=struct(Meas);


if exist('UaOutputs.m','file')==2  && ~(exist('DefineOutputs.m','file')==2)

    warning("OldInputFormat:UaOutputs","UaOutputs.m  no longer used. Rename that file to DefineOutputs.m")

end


% To reduce the size of output files, I delete some of the MUA fields
% If required this can always be recreated by the call:
%
%   MUA=UpdateMUA(CtrlVar,MUA);
%
%

MUA.DetJ=[] ; MUA.Deriv=[] ; MUA.dM=[] ;  MUA.M=[] ;  MUA.TR=[] ;    MUA.workers=[];

Nouts=nargout('DefineOutputs');
Nins=nargin('DefineOutputs');

switch Nouts

    case 0

        if Nins==6
            DefineOutputs(UserVar,CtrlVar,MUA,BCs,F,l);
        else
            UserVar=DefineOutputs(UserVar,CtrlVar,MUA,BCs,F,l,F.GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);
        end

    case 1
        if Nins==6
            DefineOutputs(UserVar,CtrlVar,MUA,BCs,F,l);
        else
            UserVar=DefineOutputs(UserVar,CtrlVar,MUA,BCs,F,l,F.GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);
        end

end

% Here the CPU and wall time tic variables are reset. The idea here is not to include CPU time used by DefineOutputs 
RunInfo.CPU.tic=duration(0,0,cputime,Format="dd:hh:mm:ss.SSS") ;
RunInfo.WallTime.tic=datetime("now",Format="dd:hh:mm:ss.SSS");


end