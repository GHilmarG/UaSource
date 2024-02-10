

%
%  Tests and examples of the use of:
%
%   LocateDetachedIslandsAndRegionsConnectedByOneNodeOnly.m
%
%

TestCase=3;  % There are currenlty 3 test cases available. 

%%  
%
%
% Note:  For TestCase=2 you will need an input file. You can get this input file from:
%
% https://livenorthumbriaac-my.sharepoint.com/:f:/g/personal/hilmar_gudmundsson_northumbria_ac_uk/EulRApgZzIlMq7mgx33DD5UB1_MDwEFIah-_CWrQcxJNnA?e=3om5cp
%
% There are several .mat files in this directory and you might want to download them all into a seperate folder and add to
% the matlab path.  These mat fiels are used in various Ua examples and tests and for plotting purposes, but none of them are required to run Ua.
%
%%

switch TestCase


    case 1
        clearvars

        CtrlVar=Ua2D_DefaultParameters();
        CtrlVar.TriNodes=3;
        CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
        CtrlVar.PlotNodalLabels=1; CtrlVar.PlotNodes=1;
        CtrlVar.PlotXYscale=1;
        CtrlVar.MeshSizeMax=0.5; CtrlVar.MeshSizeMin=0.5; CtrlVar.MeshSize=0.5;
        UserVar=[];
        CtrlVar.MeshBoundaryCoordinates=[1 NaN ; -1 -1 ; -1 -1.5 ; -1.5 -1 ; -1 -1 ; -1 0 ; 0 0.5 ; 0 0 ; -1 0 ; 1 -0.5 ; 1 -1 ; 1 -1.5 ; 0.5 -1.5 ; 1 -1 ;   ...                                                     % boundary of mesh 1
            2 NaN ; -2.0 -0.5 ; -2.0 0.5 ; -1.5 0.5 ; -1.5 -0.5 ; ...
            3 NaN ; 2.0 -0.5 ; 2.0 0.5 ; 1.5 0.5 ; 1.5 -0.5 ];

        [UserVar,MUA]=genmesh2d(UserVar,CtrlVar);

        FindOrCreateFigure("Mesh") ; PlotMuaMesh(CtrlVar,MUA); drawnow

        [Islands]=LocateDetachedIslandsAndRegionsConnectedByOneNodeOnly(CtrlVar,MUA) ;

        FindOrCreateFigure("Islands") ;
        CtrlVar.PlotNodalLabels=0;
        PlotMuaMesh(CtrlVar,MUA);
        hold on
        PlotMuaMesh(CtrlVar,MUA,Islands.Free,color="r",LineWidth=2);
        PlotMuaMesh(CtrlVar,MUA,Islands.OneNode,color="b",LineStyle="--",LineWidth=2);
        title("Islands (red) and one-node or less connection (blue)")


        %%

    case 2

        doPLots=true;

        load("EleDeactivationTestDataLarge.mat","CtrlVar","F","MUA")
        nNodesIn=MUA.Nnodes; 
        if doPLots
            FindOrCreateFigure("Mesh") ; PlotMuaMesh(CtrlVar,MUA);
        end


        fprintf("LevelSetElementDeactivation:  deactivating ALL elements downstream of calving fronts. ")
        isIceNode = F.LSF >= 0 ;              % All icy nodes
        isIceElement=AllElementsContainingGivenNodes(MUA.connectivity,isIceNode) ;  % all elements containing at least one icy node
        ElementsToBeDeactivated=~isIceElement ;

        if doPLots
            FindOrCreateFigure("MUA mesh: Ahead of level set deactivation")
            PlotMuaMesh(CtrlVar,MUA);
            hold on
            PlotMuaMesh(CtrlVar,MUA,ElementsToBeDeactivated,color="r",LineStyle="--",LineWidth=2);
            title("Level Set Deactivation")
        end

        CtrlVar.UpdateMUAafterDeactivating=false ;
        [MUA,k,l]=DeactivateMUAelements(CtrlVar,MUA,ElementsToBeDeactivated) ;

        if doPLots
            FindOrCreateFigure("MUAnew mesh: After level set deactivation")
            PlotMuaMesh(CtrlVar,MUA);
            title("After Level Set Deactivation")
        end

        CtrlVar.LocateDetachedIslandsAndRegionsConnectedByOneNodeOnly="-Islands-OneNodeOrLessConnections-" ;
        [Islands]=LocateDetachedIslandsAndRegionsConnectedByOneNodeOnly(CtrlVar,MUA) ;

        numel(find(Islands.OneNode))

        FindOrCreateFigure("Islands") ;
        CtrlVar.PlotNodalLabels=0;
        PlotMuaMesh(CtrlVar,MUA);
        hold on
        PlotMuaMesh(CtrlVar,MUA,Islands.Free,color="r",LineWidth=2);
        PlotMuaMesh(CtrlVar,MUA,Islands.OneNode,color="b",LineStyle="--",LineWidth=2);
        title("Islands (red) and one-node or less connection (blue)")

        ElementsToBeDeactivated=Islands.OneNode;
        [MUA,k,l]=DeactivateMUAelements(CtrlVar,MUA,ElementsToBeDeactivated,k,l) ;
        

        if doPLots
            FindOrCreateFigure("MUAnew mesh: After level set deactivation and island removal")
            PlotMuaMesh(CtrlVar,MUA);
            title("After Level Set Deactivation and island removal")
        end

        %%


    case 3

        % Shows how to map fields between meshes for repeated element deactivation

        CtrlVar=Ua2D_DefaultParameters();
        CtrlVar.UpdateMUAafterDeactivating=false;
        CtrlVar.TriNodes=3;
        CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
        CtrlVar.PlotNodes=1;

        CtrlVar.PlotEleLabels=1;
        CtrlVar.PlotNodalLabels=1;
        CtrlVar.PlotXYscale=1;
        CtrlVar.MeshSizeMax=0.5; CtrlVar.MeshSizeMin=0.5; CtrlVar.MeshSize=0.5;
        UserVar=[];
        CtrlVar.MeshBoundaryCoordinates=[1 NaN ; -1 -1 ; -1 -1.5 ; -1.5 -1 ; -1 -1 ; -1 0 ; 0 0.5 ; 0 0 ; -1 0 ; 1 -0.5 ; 1 -1 ; 1 -1.5 ; 0.5 -1.5 ; 1 -1 ;   ...                                                     % boundary of mesh 1
            2 NaN ; -2.0 -0.5 ; -2.0 0.5 ; -1.5 0.5 ; -1.5 -0.5 ; ...
            3 NaN ; 2.0 -0.5 ; 2.0 0.5 ; 1.5 0.5 ; 1.5 -0.5 ];

     

        [UserVar,MUA]=genmesh2d(UserVar,CtrlVar);

        nNodesIn=MUA.Nnodes ; % the original number of nodes


        FindOrCreateFigure("Mesh Original") ; PlotMuaMesh(CtrlVar,MUA);

        F.x=MUA.coordinates(:,1) ;  F.y=MUA.coordinates(:,2) ;
        UaPlots(CtrlVar,MUA,[],F.x,FigureTitle="x on original mesh") ;

        ElementsToBeDeactivated=[2 29]  ;
        [MUA,k,l]=DeactivateMUAelements(CtrlVar,MUA,ElementsToBeDeactivated)  ;

        FindOrCreateFigure("Mesh after First Deactivation") ; PlotMuaMesh(CtrlVar,MUA); 
        UaPlots(CtrlVar,MUA,[],F.x(k),FigureTitle="First deactivation") ;

        ElementsToBeDeactivated=[2 16 ]  ;  
        [MUA,k,l]=DeactivateMUAelements(CtrlVar,MUA,ElementsToBeDeactivated,k,l)  ;

        FindOrCreateFigure("Mesh after Second Deactivation") ; PlotMuaMesh(CtrlVar,MUA); 
        UaPlots(CtrlVar,MUA,[],F.x(k),FigureTitle="Second deactivation") ;


        [Islands]=LocateDetachedIslandsAndRegionsConnectedByOneNodeOnly(CtrlVar,MUA) ;

        FindOrCreateFigure("Islands") ;
        
        PlotMuaMesh(CtrlVar,MUA);
        hold on
        PlotMuaMesh(CtrlVar,MUA,Islands.Free,color="r",LineWidth=2);
        PlotMuaMesh(CtrlVar,MUA,Islands.OneNode,color="b",LineStyle="--",LineWidth=2);
        title("Islands (red) and one-node or less connection (blue)")

        ElementsToBeDeactivated=Islands.OneNode  ;  [MUA,k,l]=DeactivateMUAelements(CtrlVar,MUA,ElementsToBeDeactivated,k,l)  ;

        FindOrCreateFigure("Mesh after Third Deactivation") ; PlotMuaMesh(CtrlVar,MUA);
        UaPlots(CtrlVar,MUA,[],F.x(k),FigureTitle="Third deactivation") ;

       %  k=k1(k2(k3));                    % k gives the mapping between new and old nodal numbers, fNew=fOld(k)
       % l=1:nNodesIn ; l=l(:)+nan;       % l gives the mapping betwenen old and new, i=l(j) gives the new node number i for the old node number j
       % l(k(1:numel(k)))=1:numel(k);     %  i=l(j) is nan if the original j node was deleted


        fprintf(" done \n")




end

