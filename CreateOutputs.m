function UserVar=CreateOutputs(UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo)

warning('off','MATLAB:structOnObject')

% convert objects to structures. Otherwise it will not be possible to save/re-use these
% objects outside of Úa.
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


N=nargout('DefineOutputs');

switch N
    
    case 0
        
        DefineOutputs(UserVar,CtrlVar,MUA,BCs,F,l,F.GF);
        
    case 1
        
        UserVar=DefineOutputs(UserVar,CtrlVar,MUA,BCs,F,l,F.GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);
        
end

end