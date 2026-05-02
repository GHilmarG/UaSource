



function [UserVar,RunInfo,F1,l1,BCs1,dt]=uvh(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1)


narginchk(9,9)


[UserVar,RunInfo,F1,l1,BCs1,dt]=uvhRootFinding(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1); 

% [UserVar,RunInfo,F1,l1,BCs1,dt]=uvhLSQ(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1); 

end
