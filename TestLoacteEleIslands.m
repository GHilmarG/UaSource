

%%

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