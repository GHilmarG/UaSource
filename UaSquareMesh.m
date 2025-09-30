

function [coordinates,connectivity]=UaSquareMesh(CtrlVar)

%%

if ~isfield(CtrlVar,"UaSquareMesh") ...
        || ~isfield(CtrlVar.UaSquareMesh,"xmin")  ...
        || ~isfield(CtrlVar.UaSquareMesh,"xmax")  ...
        ||~isfield(CtrlVar.UaSquareMesh,"nx")

    fprintf("When using the UaSquareMesh mesh generator, the fields : \n ")
    fprintf(" CtrlVar.UaSquareMesh.xmin \n CtrlVar.UaSquareMesh.xmax \n CtrlVar.UaSquareMesh.ymin \n CtrlVar.UaSquareMesh.ymax \n CtrlVar.UaSquareMesh.nx \n CtrlVar.UaSquareMesh.ny \n")
    fprintf("Must all be defined. \n ")
    
  
    error("UaSquareMesh:Inputs","not all input fields defined")
  

end

xmin=CtrlVar.UaSquareMesh.xmin;
xmax=CtrlVar.UaSquareMesh.xmax;
ymin=CtrlVar.UaSquareMesh.ymin;
ymax=CtrlVar.UaSquareMesh.ymax;

nx=CtrlVar.UaSquareMesh.nx;
ny=CtrlVar.UaSquareMesh.ny;


% xmin=-10 ; xmax=10 ; ymin=-5 ; ymax=5;  dx=1 ; dy=1 ;  nx=round((xmax-xmin)/dx); ny=round((ymax-ymin)/dy); 


x=linspace(xmin,xmax,nx+1);
y=linspace(ymin,ymax,ny+1);
[X,Y]=ndgrid(x,y);


x=X(:); 
y=Y(:) ;
DT = delaunayTriangulation(x,y) ; 

connectivity=DT.ConnectivityList;
coordinates=DT.Points ; 

% figure(999) ; plot(x,y,".") ; axis equal
% figure(1002); triplot(DT) ; axis equal
% CtrlVar=Ua2D_DefaultParameters(); MUA=CreateMUA(CtrlVar,DT.ConnectivityList,DT.Points);

%%