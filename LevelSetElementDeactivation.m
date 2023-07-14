function ElementsToBeDeactivated=LevelSetElementDeactivation(RunInfo,CtrlVar,MUA,F,ElementsToBeDeactivated)


if isempty(F.LSF) ; return ; end




if CtrlVar.LevelSetMethodAutomaticallyDeactivateElements

    if CtrlVar.LevelSetEvolution=="-Prescribed-"

        % get rid of ALL elements downstream of the calving fronts


        fprintf("LevelSetElementDeactivation:  deactivating ALL elements downstream of calving fronts. ")
        isIceNode = F.LSF >= 0 ;              % All icy nodes
        isIceElement=AllElementsContainingGivenNodes(MUA.connectivity,isIceNode) ;  % all elements containing at least one icy node
        ElementsToBeDeactivated=~isIceElement ;

    else


        CtrlVar.LineUpGLs=false ;

        if isnan(CtrlVar.LevelSetMethodStripWidth)

            fprintf("LevelSetElementDeactivation: The variable CtrlVar.LevelSetMethodStripWidth is not defined.\n")
            fprintf("LevelSetElementDeactivation: Elements will not be deactivated based on the value of the level set.\n")
            warning("LevelSetEquation:ParameterNotDefined","The variable CtrlVar.LevelSetMethodStripWidth needs to be defined.")

        else

            fprintf("LevelSetElementDeactivation:  deactivating ALL elements downstream of calving fronts. ")


            % I only deactivate elements if:
            % 1) all nodes are have LSF value below threhold
            % 2) all nodes have a negative LSF value

            %

            [xc,yc]=CalcMuaFieldsContourLine(CtrlVar,MUA,F.LSF,0);
            DistNod=pdist2([xc(:) yc(:)],[F.x F.y],'euclidean','Smallest',1) ;

            DistNod=DistNod(:) ;

            isIcyNode = DistNod < CtrlVar.LevelSetMethodStripWidth  | F.LSF > 0 ;
            isIceElement=AllElementsContainingGivenNodes(MUA.connectivity,isIcyNode) ;
            ElementsToBeDeactivated=~isIceElement ;


        end

    end
end

% UaPlots(CtrlVar,MUA,F,F.h) ; hold on ;  PlotMuaMesh(CtrlVar,MUA) ; hold on ; PlotMuaMesh(CtrlVar,MUA,ElementsToBeDeactivated,"r") ;

