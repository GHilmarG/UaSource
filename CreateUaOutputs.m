function UserVar=CreateUaOutputs(UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo)

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


N=nargout('UaOutputs');

switch N
    
    case 0
        
        UaOutputs(UserVar,CtrlVar,MUA,BCs,UaVars,l,GF);
        
    case 1
        
        UserVar=UaOutputs(UserVar,CtrlVar,MUA,BCs,F,l,F.GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo);
        
end

end