function [UserVar,RunInfo,LSF,l,LSFqx,LSFqy,BCs]=LevelSetEquationAnalyticalInitialisation(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l)

nargoutchk(7,7)


LSFqx=[] ; LSFqy=[] ; Mask=[] ; 

% Don't redefine F0.LSF as F1.LSF, doing so would push the solution back in tiem
F1.LSF=F0.LSF ;
Threshold=0 ;    % Level Set value

% Here F0.LSF is the original, and F1.LSF will be the re-initilized LSF
% fix the LSF field for all nodes of elements around the level.

if CtrlVar.LevelSetInitBCsZeroLevel
    % Use BCs to fix the level set over all elements that the level
    % goes through. This ensures that the level can not shift during
    % initialisation.
    Mask=CalcMeshMask(CtrlVar,MUA,F0.LSF,Threshold);
    
    LSFFixedNodeUnmodified=BCs.LSFFixedNode ;
    LSFFixedValueUnmodified=BCs.LSFFixedValue ;
    
    BCs.LSFFixedNode= [LSFFixedNodeUnmodified ; find(Mask.NodesOn)];  % add the nodes of the "On" elements, ie all elements containing the zero level
    BCs.LSFFixedValue=[LSFFixedValueUnmodified ; F0.LSF(Mask.NodesOn) ];
   % 

end

% CtrlVar.LevelSetReinitializePDist=1;
if  CtrlVar.LevelSetReinitializePDist

    %% After having located the 0 level, now do a rough re-initialisation using signed distance function. After this I then do a full
    % non-linear FAB solve with the level-set fixed as boundary conditions on the LSF.
    % This will in most cases not be needed, but

    if  isfield(CtrlVar,'CtrlVar.LevelSetTestString') &&  contains(CtrlVar.LevelSetTestString,"-xc/yc nodes-")
        xC=F0.x(Mask.NodesOn) ; yC=F0.y(Mask.NodesOn) ;
    else
        CtrlVar.LineUpGLs=false ;
        [xC,yC]=CalcMuaFieldsContourLine(CtrlVar,MUA,F0.LSF,Threshold);
    end

    % It should be OK to do this with LSF at both 0 and 1 as I have already
    % found the location of the level set for F0 and this will be enforced
    % throught the BCs.
    [LSF,UserVar,RunInfo]=SignedDistUpdate(UserVar,RunInfo,CtrlVar,MUA,F0.LSF,xC,yC);
    F0.LSF=LSF ;
    F1.LSF=LSF ;
end
%%


[UserVar,RunInfo,LSF,l]=LevelSetFixPointSolver(UserVar,RunInfo,CtrlVar,MUA,BCs,F0,F1,l) ; 
% var her
return


end