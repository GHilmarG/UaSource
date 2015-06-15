%   1     2      3        4        5      6     7     8      9
X=[0 0 ; 1 0 ; 1 0.25 ; 1 0.5 ; 1 0.75 ; 1 1 ; 0 1 ; 2 0 ; 2 1];
edge=[1 2 ; 2 3 ; 3 4 ; 4 5 ; 5 6 ; 6 7 ; 7 1 ; 2  8 ; 8 9 ; 9 6 ];


hdata.MeshSizeMax=1;
options.output=true;
%[p,t,stats] = mesh2d(X,edge) ; %,hdata,options);
face{1}=[1;2;3;4;5;6;7];
[coordinates,connectivity,fnum,stats]=meshfaces(X,edge,[],hdata,options);

alphaSweep=0.01;

[coordinates,connectivity] = ElementSweep(coordinates,connectivity,alphaSweep);
[coordinates,connectivity] = NodalSweep(coordinates,connectivity,alphaSweep);


figure ; triplot(connectivity,coordinates(:,1),coordinates(:,2))

CtrlVar.PlotLabels=1;
PlotFEmesh(coordinates,connectivity,CtrlVar)
    

%%%


%   1     2      3        4        5      6     7     8      9
X=[0 0 ; 1 0 ; 1 0.25 ; 1 0.5 ; 1 0.75 ; 1 1 ; 0 1 ; 2 0 ; 2 1];

%      1     2     3      4    5     6     7     8      9      10   11
edge=[1 2 ; 2 3 ; 3 4 ; 4 5 ; 5 6 ; 6 7 ; 7 1 ; 2  8 ; 8 9 ; 9 6 ; 6 2];



hdata.MeshSizeMax=1;
options.output=true;
%[p,t,stats] = mesh2d(X,edge) ; %,hdata,options);
face{1}=[1;2;3;4;5;6;7];
face{2}=[8;9;10;11];
[coordinates,connectivity,fnum,stats]=meshfaces(X,edge,face,hdata,options);

alphaSweep=0.01;

[coordinates,connectivity] = ElementSweep(coordinates,connectivity,alphaSweep);
[coordinates,connectivity] = NodalSweep(coordinates,connectivity,alphaSweep);


figure ; triplot(connectivity,coordinates(:,1),coordinates(:,2))

CtrlVar.PlotLabels=1;
PlotFEmesh(coordinates,connectivity,CtrlVar)
    

%%%



