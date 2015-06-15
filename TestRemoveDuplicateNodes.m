

%%
CtrlVar.PlotNodes=1; CtrlVar.PlotLabels=1;

%
coordinates=[0 0 ; 1 0 ; 2 0 ; 0 1 ; 1 1];
connectivity=[1 4 2 ; 2 4 5 ; 2 5 3];

figure
PlotFEmesh(coordinates,connectivity,CtrlVar)

[coordinates,connectivity]=RemoveDuplicateNodes(coordinates,connectivity);



% duplicate element
coordinates=[0 0 ; 1 0 ; 2 0 ; 0 1 ; 1 1];
connectivity=[1 4 2 ; 2 4 5 ; 2 5 3 ; 3 2 5];

figure
PlotFEmesh(coordinates,connectivity,CtrlVar)

[coordinates,connectivity]=RemoveDuplicateNodes(coordinates,connectivity);
coordinates
connectivity




%% duplicate node and elements
coordinates=[0 0 ; 1 0 ; 2 0 ; 0 1 ; 1 1 ; 1 1 ];
connectivity=[1 4 2 ; 2 4 5 ; 2 5 3 ; 3 2 5];

figure ; PlotFEmesh(coordinates,connectivity,CtrlVar)

[coordinates,connectivity]=RemoveDuplicateNodes(coordinates,connectivity);
figure ; PlotFEmesh(coordinates,connectivity,CtrlVar)
coordinates
connectivity

%%
%% duplicate node and elements
coordinates=[0 0 ; 1 0 ; 2 0 ; 0 1 ; 1 1 ; 1 1 ; 0.5 0.5];
connectivity=[1 4 2 ; 2 4 5 ; 2 5 3 ; 3 2 5];

figure
PlotFEmesh(coordinates,connectivity,CtrlVar)

[coordinates,connectivity]=RemoveDuplicateNodes(coordinates,connectivity);
figure ; PlotFEmesh(coordinates,connectivity,CtrlVar)
coordinates
connectivity




