function [UserVar,RunInfo,F1,l1,BCs1,dt]=uvh(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1)

nargoutchk(6,6)

% At the moment, just a wrapper



[UserVar,RunInfo,F1,l1,BCs1,dt]=FIuvh2D(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1);



end