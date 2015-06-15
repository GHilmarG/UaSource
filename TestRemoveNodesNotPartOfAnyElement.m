%% TestRemoveNodesNotPartOfAnyElement(coordinates,connectivity)


close all

coordinates=[1.5 0  ; 0 0 ; 0 0.5 ; 1 0 ; 2 0 ; 0 1 ; 1 1 ; 1 1 ; 0.5 0.5 ; 0 0];
connectivity=[8 5 4 ; 2 6 4 ; 4 6 7 ; 4 8 5 ];

%%
CtrlVar.PlotNodes=1; CtrlVar.PlotLabels=1; 
figure ; PlotFEmesh(coordinates,connectivity,CtrlVar)
%%

[coordinates,connectivity]=RemoveNodesNotPartOfAnyElement(coordinates,connectivity);


coordinates
connectivity


figure ; PlotFEmesh(coordinates,connectivity,CtrlVar)

[coordinates,connectivity]=RemoveDuplicateNodes(coordinates,connectivity);

figure ; PlotFEmesh(coordinates,connectivity,CtrlVar)

coordinates
connectivity




