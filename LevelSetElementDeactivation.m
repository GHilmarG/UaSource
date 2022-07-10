function ElementsToBeDeactivated=LevelSetElementDeactivation(RunInfo,CtrlVar,MUA,F,ElementsToBeDeactivated)


if isempty(F.LSF) ; return ; end





% Mask=CalcMeshMask(CtrlVar,MUA,F.LSF,0);
%      CtrlVar.LevelSetMethodAutomaticallyDeactivateElementsThreshold=CtrlVar.LevelSetMethodStripWidth/2 ;
%      Inode=F.LSF<CtrlVar.LevelSetMethodAutomaticallyDeactivateElementsThreshold ;
%      Iele=MuaElementsContainingGivenNodes(CtrlVar,MUA,find(Inode),Mask.ElementsOut,"all") ;
%      ElementsToBeDeactivated=ElementsToBeDeactivated | Iele ;



if CtrlVar.LevelSetMethodSolveOnAStrip



    CtrlVar.LineUpGLs=false ; Threshold=0 ;

    [xc,yc]=CalcMuaFieldsContourLine(CtrlVar,MUA,F.LSF,Threshold);
    DistNod=pdist2([xc(:) yc(:)],[F.x F.y],'euclidean','Smallest',1) ;

    DistNod=DistNod(:) ;
    DistEle=Nodes2EleMean(MUA.connectivity,DistNod) ;  % note, this is now an element-valued distance function
    LSFEle=Nodes2EleMean(MUA.connectivity,F.LSF) ;

    if isnan(CtrlVar.LevelSetMethodStripWidth)

        fprintf("The variable CtrlVar.LevelSetMethodStripWidth needs to be defined.\n")
        error("LevelSetEquation:ParameterNotDefined","The variable CtrlVar.LevelSetMethodStripWidth needs to be defined.")

    end


    %     NegLSFNode=F.LSF < 0;
    %     OutsideStripNodes=DistNode<CtrlVar.LevelSetMethodStripWidth ;
    %     Inode=NegLSFNode & OutsideStripNodes ;
    %     Iele=MuaElementsContainingGivenNodes(CtrlVar,MUA,Inode,[],"all") ;

    ElementsToBeDeactivated=DistEle>CtrlVar.LevelSetMethodStripWidth & LSFEle < 0 ; % distance downstream of calving front


end


% figure ; UaPlots(CtrlVar,MUA,F,F.h) ; hold on ;  PlotMuaMesh(CtrlVar,MUA) ; hold on ; PlotMuaMesh(CtrlVar,MUA,ElementsToBeDeactivated,"r") ;

