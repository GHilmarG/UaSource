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

CtrlVar.PlotLabels=1;CtrlVar.PlotNodes=1 ; 
PlotFEmesh(coordinates,connectivity,CtrlVar)
    

%%


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
    



%%


%   1     2      3    4       5      6       7       8        9     10           11
X=[0 0 ; 2 0 ; 2 1 ; 0 1 ; 1 0. ; 1 0.3 ; 1 0.4 ; 1 0.5 ; 1 0.6 ; 0.9 0.6 ; 0.8 0.6 ];

%      1     2     3     4    5     6       
edge=[1 2 ; 2 3 ; 3 4 ; 4 1 ; 5 6 ; 6 7 ; 7 8 ; 8 9; 9 10  ; 10  11 ; 11 10 ;  10  9 ; 9 8 ; 8 7 ; 7 6 ; 6 5 ]; 


hdata.MeshSizeMax=1;
options.output=true;
%[p,t,stats] = mesh2d(X,edge) ; %,hdata,options);
[coordinates,connectivity,fnum,stats]=meshfaces(X,edge,[],[],options);

alphaSweep=0.01;

[coordinates,connectivity] = ElementSweep(coordinates,connectivity,alphaSweep);
[coordinates,connectivity] = NodalSweep(coordinates,connectivity,alphaSweep);


figure ; triplot(connectivity,coordinates(:,1),coordinates(:,2))

CtrlVar.PlotLabels=1;
PlotFEmesh(coordinates,connectivity,CtrlVar)
    

%%

%%


%   1     2      3    4       5      6       7       8        9     10           11
X=[0 0 ; 2 0 ; 2 1 ; 0 1 ; 1 -1. ; 1 0.3 ; 1 0.4 ; 1 0.5 ; 1 0.6 ; 0.9 0.6 ; 0.8 1 ; 1.3 0.25; 1.6 0.25 ; 1.6 0.75 ;1.3 0.75 ];

%      1     2     3     4    5     6       
edge=[1 2 ; 2 3 ; 3 4 ; 4 1 ; 5 6 ; 6 7 ; 7 8 ; 8 9; 9 10  ; 10  11 ; 11 10 ;  10  9 ; 9 8 ; 8 7 ; 7 6 ; 6 5 ; 12 13 ; 13 14 ; 14 15 ; 15 12]; 


hdata.MeshSizeMax=1;
options.output=true;
%[p,t,stats] = mesh2d(X,edge) ; %,hdata,options);
[coordinates,connectivity,fnum,stats]=meshfaces(X,edge,[],[],options);

alphaSweep=0.01;

[coordinates,connectivity] = ElementSweep(coordinates,connectivity,alphaSweep);
[coordinates,connectivity] = NodalSweep(coordinates,connectivity,alphaSweep);


figure ; triplot(connectivity,coordinates(:,1),coordinates(:,2))

CtrlVar.PlotLabels=1;
PlotFEmesh(coordinates,connectivity,CtrlVar)
    


%%


%   1     2      3    4       5      6       7       8        9     10           11



y=0.01:0.01:0.99; x=1+0.1*sin(2*y*2*pi);
X=[0 0 ; 2 0 ; 2 1 ; 0 1 ; [x' y']];

ea=[4+(1:(numel(x)-1)), 4+numel(x):-1:6 ]
eb=CircShiftVector(ea,1);


edge=[1 2 ; 2 3 ; 3 4 ; 4 1 ; [ea' eb'] ];


hdata.MeshSizeMax=1;
options.output=true;
%[p,t,stats] = mesh2d(X,edge) ; %,hdata,options);
[coordinates,connectivity,fnum,stats]=meshfaces(X,edge,[],hdata,options);


alphaSweep=0.01;

[coordinates,connectivity] = ElementSweep(coordinates,connectivity,alphaSweep);
[coordinates,connectivity] = NodalSweep(coordinates,connectivity,alphaSweep);


figure ; triplot(connectivity,coordinates(:,1),coordinates(:,2))

CtrlVar.PlotLabels=0;
PlotFEmesh(coordinates,connectivity,CtrlVar)
    


%%
% moving the mesh
%   1     2      3    4       5      6       7       8        9     10           11
close all ; clear all 
y=.1:0.01:0.9; x=1+0.1*sin(2*y*2*pi);
X=[0 0 ; 2 0 ; 2 1 ; 0 1 ; [x' y']];

ea=[4+(1:(numel(x)-1)), 4+numel(x):-1:6 ];
eb=CircShiftVector(ea,1);
edge=[1 2 ; 2 3 ; 3 4 ; 4 1 ; [ea' eb'] ];

hdata.MeshSizeMax=1;
options.output=true;
%[p,t,stats] = mesh2d(X,edge) ; %,hdata,options);
[coordinates,connectivity,fnum,stats]=meshfaces(X,edge,[],hdata,options);
alphaSweep=0.01;
[coordinates,connectivity] = ElementSweep(coordinates,connectivity,alphaSweep);
[coordinates,connectivity] = NodalSweep(coordinates,connectivity,alphaSweep);
figure ; triplot(connectivity,coordinates(:,1),coordinates(:,2))

CtrlVar.PlotLabels=0; CtrlVar.PlotNodes=1; 
PlotFEmesh(coordinates,connectivity,CtrlVar)
    

figure
ea=[4+(1:(numel(x)))];

Constraints=[1 2 ; 2 3 ; 3 4 ; 4 1 ; [ea(1:end-1)' ea(2:end)'] ];
dt = DelaunayTri(X(:,1),X(:,2),Constraints);
clf; triplot(dt); axis equal;
descriptors.tri = dt.pointLocation(coordinates(:,1), coordinates(:,2));
descriptors.baryCoords = dt.cartToBary(descriptors.tri, coordinates);

% GL moves

x=1.1+0.1*sin(2*y*2*pi);
X=[0 0 ; 2 0 ; 2 1 ; 0 1 ; [x' y']];
tr = TriRep(dt(:,:),X(:,1),X(:,2));
figure
triplot(tr); axis equal tight
title(' GL moves')


Xnew = tr.baryToCart(descriptors.tri, descriptors.baryCoords);
tr = TriRep(connectivity, Xnew);
figure; triplot(tr);
axis equal tight;
xlabel('Morphed Mesh', 'fontweight','b');


%%
% moving the mesh (two faces)
%   1     2      3    4       5      6       7       8        9     10           11
close all ; clear all 
y=0:0.02:1; x=1+0.1*sin(2*y*2*pi);
X=[0 0 ; 2 0 ; 2 1 ; 0 1 ; [x' y']];
iGL=[4+(1:(numel(x)))];
GLedge=[iGL(1:end-1)' iGL(2:end)'];

edge=[1 iGL(1) ; iGL(1) 2 ; 2 3 ; 3 iGL(end) ; iGL(end) 4 ; 4 1 ; GLedge ];
face{1}=[1;[7:6+size(GLedge,1)]';5;6];
face{2}=[2;3;4;[6+size(GLedge,1):-1:7]'];


hdata.MeshSizeMax=0.1;
options.output=true;
%[p,t,stats] = mesh2d(X,edge) ; %,hdata,options);
[coordinates,connectivity,fnum,stats]=meshfaces(X,edge,face,hdata,options);


alphaSweep=0.01;
[coordinates,connectivity] = ElementSweep(coordinates,connectivity,alphaSweep);
[coordinates,connectivity] = NodalSweep(coordinates,connectivity,alphaSweep);
figure ; triplot(connectivity,coordinates(:,1),coordinates(:,2))

CtrlVar.PlotLabels=0; CtrlVar.PlotNodes=1; 
PlotFEmesh(coordinates,connectivity,CtrlVar)
    

figure


Constraints=[1 2 ; 2 3 ; 3 4 ; 4 1 ; GLedge];
dt = DelaunayTri(X(:,1),X(:,2),Constraints);
clf; triplot(dt); axis equal;
descriptors.tri = dt.pointLocation(coordinates(:,1), coordinates(:,2)); 
descriptors.baryCoords = dt.cartToBary(descriptors.tri, coordinates);

%% GL moves
%
%close all

%%

doVideo=0;

close all

if doVideo==1
    scrsz = get(0,'ScreenSize');
    figure('Position',[50 50 scrsz(3)/1.1 scrsz(4)/1.1])
else
    figure
end



if doVideo==1
    vidObj = VideoWriter('MeshMorphing.avi');
    open(vidObj);
    
end

set(gca,'nextplot','replacechildren');

N=5;
for I=1:N
    I
    %subplot(N,1,I)
    
    
    
    x=I/N/25+1+0.1*sin(2*y*2*pi);
    X=[0 0 ; 2 0 ; 2 1 ; 0 1 ; [x' y']];
    
    %     figure
    %     triplot(tr); axis equal tight
    %     title(' GL moves')
    
    tr = TriRep(dt(:,:),X(:,1),X(:,2));
    coordinates = tr.baryToCart(descriptors.tri, descriptors.baryCoords);
    
    %[coordinates,connectivity6]=tri3to6(coordinates,connectivity);
    %PlotFEmesh(coordinates,connectivity6,CtrlVar);
    tr = TriRep(connectivity, coordinates);
    triplot(tr)
    hold on
    plot(x,y,'r','LineWidth',2)
    axis equal tight
    if doVideo==1
        currFrame = getframe;
        writeVideo(vidObj,currFrame);
    end
    hold off
end

if doVideo==1
    close(vidObj)
end


