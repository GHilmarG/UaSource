function ElementsToBeDeactivated=LevelSetElementDeactivation(RunInfo,CtrlVar,MUA,F,ElementsToBeDeactivated)
         
    
    if isempty(F.LSF) ; return ; end
    
    
    
    
    
    Mask=CalcMeshMask(CtrlVar,MUA,F.LSF,0);
    
    % Iele=Mask.ElementsOut; 

    % This approach gives one extra row of elements downstream of calving fronts
     CtrlVar.LevelSetMethodAutomaticallyDeactivateElementsThreshold=CtrlVar.LevelSetMethodStripWidth/2 ;
     Inode=F.LSF<CtrlVar.LevelSetMethodAutomaticallyDeactivateElementsThreshold ;
     Iele=MuaElementsContainingGivenNodes(CtrlVar,MUA,find(Inode),Mask.ElementsOut,"all") ;


    ElementsToBeDeactivated=ElementsToBeDeactivated | Iele ; 
    
end
