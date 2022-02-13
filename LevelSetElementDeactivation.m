function ElementsToBeDeactivated=LevelSetElementDeactivation(RunInfo,CtrlVar,MUA,F,ElementsToBeDeactivated)
         
    
    if isempty(F.LSF) ; return ; end
    
    
    
    
    
    Mask=CalcMeshMask(CtrlVar,MUA,F.LSF,0);
    
    % Iele=Mask.ElementsOut; 

   
%      CtrlVar.LevelSetMethodAutomaticallyDeactivateElementsThreshold=CtrlVar.LevelSetMethodStripWidth/2 ;
%      Inode=F.LSF<CtrlVar.LevelSetMethodAutomaticallyDeactivateElementsThreshold ;
%      Iele=MuaElementsContainingGivenNodes(CtrlVar,MUA,find(Inode),Mask.ElementsOut,"all") ;
%      ElementsToBeDeactivated=ElementsToBeDeactivated | Iele ; 


     ElementsToBeDeactivated=Mask.ElementsOut; 

    
end
