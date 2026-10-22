



function [UserVar,RunInfo,F1,l1,BCs1,dt]=uvh(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1)


narginchk(9,9)


if ~isfield(CtrlVar,"CminWidth") || ~isfield(CtrlVar,"AGlenminWidth")
    CtrlVar=SetSmoothFloorWidths(CtrlVar,F0);
end


[UserVar,RunInfo,F1,l1,BCs1,dt]=uvhRootFinding(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1);



return

end
