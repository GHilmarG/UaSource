function [coordinates,connectivity,MUA]=UaSquareMesh(CtrlVar)

%% Simple mesh generator for square domains.
%
% The idea is that the mesh is completely regular. This can not be achieved with gmsh or mesh2d.
%
% This is good for, for example, convergence studies, and when applying periodic boundary conditions, which is often done
% using very simple mesh geometries, e.g. square shaped domains. 
%
%
% To arrive at a mesh consisting of regular smaller squares divided into four triangles, set
%
%   CtrlVar.UaSquareMesh.Refine=true;
% 
% Otherwise the mesh will be based on straightforward delaunay triangulation which, typically, results in a mesh that is not
% symmetrical with respect to a rotation about 90 degrees, and mirror reflections around x and y.
%
%
% This simple mesh generator is selected by setting 
%
%
%    CtrlVar.MeshGenerator="UaSquareMesh"; 
%
% in DefineInitialInputs.m
%
% The extent of the square is specified by setting
%
%   CtrlVar.UaSquareMesh.xmin  
%   CtrlVar.UaSquareMesh.xmax 
%   CtrlVar.UaSquareMesh.ymin   
%   CtrlVar.UaSquareMesh.ymax
%
%
% And this square is then subdivided into nx and ny elements in the x and the y directions.
%
%   CtrlVar.UaSquareMesh.nx;
%   CtrlVar.UaSquareMesh.ny;
%
% If additionally, and this is recommended, 
%
%   CtrlVar.UaSquareMesh.Refine=true;
%
% each triangle is additionally subdivided from the longest vertex.
%
% If, for example, 
%
%
%   CtrlVar.UaSquareMesh.xmin=-50e3 ;
%   CtrlVar.UaSquareMesh.xmax=50e3 ;
%   CtrlVar.UaSquareMesh.ymin=-50e3;
%   CtrlVar.UaSquareMesh.ymax=50e3;
% 
%   CtrlVar.UaSquareMesh.nx=10;
%   CtrlVar.UaSquareMesh.ny=10;
%
% and
%
% CtrlVar.UaSquareMesh.Refine=false;
%
% then the resulting size of every element will be:
%
% 100e3/10=10
%
% If 
%
%   CtrlVar.UaSquareMesh.Refine=true;
% 
% then the element size will be
%
% 100e3/10/sqrt(2) =   7071.06781186548
%
%%

if ~isfield(CtrlVar,"UaSquareMesh") ...
        || ~isfield(CtrlVar.UaSquareMesh,"xmin")  ...
        || ~isfield(CtrlVar.UaSquareMesh,"xmax")  ...
        ||~isfield(CtrlVar.UaSquareMesh,"nx") ...
        ||~isfield(CtrlVar.UaSquareMesh,"nx")

    fprintf("When using the UaSquareMesh mesh generator, the fields : \n ")
    fprintf(" CtrlVar.UaSquareMesh.xmin \n CtrlVar.UaSquareMesh.xmax \n CtrlVar.UaSquareMesh.ymin \n CtrlVar.UaSquareMesh.ymax \n CtrlVar.UaSquareMesh.nx \n CtrlVar.UaSquareMesh.ny \n")
    fprintf("Must all be defined. \n ")


    error("UaSquareMesh:Inputs","not all input fields defined")


end


if ~isfield(CtrlVar.UaSquareMesh,"Refine")
    CtrlVar.UaSquareMesh.Refine=true;
end

if ~isfield(CtrlVar,"TriNodes")
    CtrlVar.TriNodes=3;
end

if ~isfield(CtrlVar,"QueadRules2021")
    CtrlVar.QuadRules2021=true;
end



xmin=CtrlVar.UaSquareMesh.xmin;
xmax=CtrlVar.UaSquareMesh.xmax;
ymin=CtrlVar.UaSquareMesh.ymin;
ymax=CtrlVar.UaSquareMesh.ymax;

nx=CtrlVar.UaSquareMesh.nx;
ny=CtrlVar.UaSquareMesh.ny;

% If the user only defines nx, than calculate a reasonable ny
if isfinite(nx) && isnan(ny)
    ny=round((CtrlVar.UaSquareMesh.ymax-CtrlVar.UaSquareMesh.ymin)/(CtrlVar.UaSquareMesh.xmax-CtrlVar.UaSquareMesh.xmin))*CtrlVar.UaSquareMesh.nx;
end

% If the user only defines nx, than calculate a reasonable ny
if isfinite(ny) && isnan(nx)
    nx=round((CtrlVar.UaSquareMesh.xmax-CtrlVar.UaSquareMesh.xmin)/(CtrlVar.UaSquareMesh.ymax-CtrlVar.UaSquareMesh.ymin))*CtrlVar.UaSquareMesh.ny;
end

% xmin=-10 ; xmax=10 ; ymin=-5 ; ymax=5;  dx=1 ; dy=1 ;  nx=round((xmax-xmin)/dx); ny=round((ymax-ymin)/dy);


x=linspace(xmin,xmax,nx+1);
y=linspace(ymin,ymax,ny+1);
[X,Y]=ndgrid(x,y);


x=X(:);
y=Y(:) ;
DT = delaunayTriangulation(x,y) ;

connectivity=DT.ConnectivityList;
coordinates=DT.Points ;


if CtrlVar.UaSquareMesh.Refine

    % For the refinement, I need to temporarily create a MUA structure. 
    CtrlVar.CalcMUA_Derivatives=false;
    CtrlVar.FindMUA_Boundary=false;
    CtrlVar.MUA.MassMatrix=false ;
    CtrlVar.MUA.StiffnessMatrix=false;
    CtrlVar.MUA.DecomposeMassMatrix=false ;
    CtrlVar.MUA.DecomposeMassMatrix=false ;
    CtrlVar.Parallel.uvAssembly.spmd.isOn=false ;
    CtrlVar.Parallel.uvhAssembly.spmd.isOn=false ;
    CtrlVar.InfoLevelAdaptiveMeshing=0 ;

    if ~isfield(CtrlVar,"AdaptMeshUntilChangeInNumberOfElementsLessThan")
        CtrlVar.AdaptMeshUntilChangeInNumberOfElementsLessThan=0;
    end

    if ~isfield(CtrlVar,"MaxNumberOfElements")
        CtrlVar.MaxNumberOfElements=inf;
    end

    MUA=CreateMUA(CtrlVar,connectivity,coordinates);


    if  isfield(CtrlVar,"MeshBoundaryCoordinates") &&  ~isempty(CtrlVar.MeshBoundaryCoordinates)

        % if MeshBoundaryCoordinates have been defined, eliminate elements outside of the desired computational boundary
        xy=[MUA.xEle,MUA.yEle];
        isInside=InsideOutside(xy,CtrlVar.MeshBoundaryCoordinates) ;
        CtrlVar.UpdateMUAafterDeactivating=true;
        MUA=DeactivateMUAelements(CtrlVar,MUA,~isInside) ;

    end


    ElementsToBeRefined=true(MUA.Nele,1);             % refine all elements
    ElementsToBeCoarsened=false(MUA.Nele,1);          % unrefine none
    RunInfo=[];

    CtrlVar.MeshRefinementMethod="explicit:local:newest vertex bisection";
    %CtrlVar.MeshRefinementMethod="explicit:local:red-green"; CtrlVar.LocalAdaptMeshSmoothingIterations=0;

    MUA=LocalMeshRefinement(CtrlVar,RunInfo,MUA,ElementsToBeRefined,ElementsToBeCoarsened);

    coordinates=MUA.coordinates;
    connectivity=MUA.connectivity ;

end



[Nele,Nod]=size(connectivity);

if Nele > CtrlVar.MaxNumberOfElements
    fprintf("UaSquareMesh: Too many elements! \n")
    fprintf("  The maximum allowed number of elements is CtrlVar.MaxNumberOfElements=%i \n",CtrlVar.MaxNumberOfElements)
    fprintf("  This is larger than the number of elments in the mesh which is %i \n",Nele)
    error("UaSquareMesh:TooManyElements","Too Many Elements")
end


end

% figure(999) ; plot(x,y,".") ; axis equal
% figure(1002); triplot(DT) ; axis equal
% CtrlVar=Ua2D_DefaultParameters(); MUA=CreateMUA(CtrlVar,DT.ConnectivityList,DT.Points);

%%