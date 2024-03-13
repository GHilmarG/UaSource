function UserVar=CreateOutputs(UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo)

warning('off','MATLAB:structOnObject')

% convert objects to structures. Otherwise it will not be possible to save/re-use these
% objects outside of �a.
%
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

end