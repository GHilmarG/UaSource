function [UserVar,MUA]=genmesh2d(UserVar,CtrlVar,MeshBoundaryCoordinates,edge,face,GmshBackgroundScalarField)

%% Generate FE mesh
%
% [MUA,FEmeshTriRep]=genmesh2d(CtrlVar,MeshBoundaryCoordinates,edge,face,GmshBackgroundScalarField)
%
% FEmeshTriRep is an instance of matlab TriRep class. Only based on the corner nodes of
%    triangles! Interior and edge nodes of higher order elements missing

nargoutchk(2,2)


options.output=false;

switch lower(CtrlVar.MeshGenerator)
    
    case 'mesh2d'
        hdata.fun = @hfun;
        hdata.args  = {CtrlVar.MeshSize};
        
        if nargin<3 || isempty(edge) || isempty(face)
            [coordinates,connectivity] = mesh2d(MeshBoundaryCoordinates,[],hdata,options);
        else
            [coordinates,connectivity]=meshfaces(MeshBoundaryCoordinates,edge,face,hdata,options);
        end
        
    case 'gmsh'
        
        if nargin<5
            GmshBackgroundScalarField=[];
        end

        [coordinates,connectivity]=GmshInterfaceRoutine(CtrlVar,MeshBoundaryCoordinates,GmshBackgroundScalarField);
        
    otherwise
        error('Mesh generator not correctly defined. Define variable CtrlVar.MeshGenerator {mesh2d|gmsh} ')
end


[coordinates,connectivity]=ChangeElementType(coordinates,connectivity,CtrlVar.TriNodes);

if CtrlVar.sweep
    [coordinates,connectivity] = ElementSweep(coordinates,connectivity,CtrlVar.SweepAngle);
    [coordinates,connectivity] = NodalSweep(coordinates,connectivity,CtrlVar.SweepAngle);
end


if CtrlVar.CuthillMcKee
    M=connectivity2adjacency(connectivity);
    [coordinates,connectivity] = CuthillMcKeeFE(coordinates,connectivity,M);
end


% removing possible duplicate nodes
% There should in principle be no need to do this
% but this is fast so not much time is lost.
% tolerance=100*eps;
% [coordinates,connectivity]=RemoveDuplicateNodes(coordinates,connectivity,tolerance);

connectivity=TestAndCorrectForInsideOutElements(CtrlVar,coordinates,connectivity);

%% Possible user modifications to coordinates and connectivity
[UserVar,coordinates,connectivity]=GetMeshModifications(UserVar,CtrlVar,coordinates,connectivity);

% [K,~,J]=unique(connectivity(:));
% connectivity=reshape(J,size(connectivity));
% coordinates=coordinates(K,:);

%FEmeshTriangulation=CreateFEmeshTriRep(connectivity,coordinates);

MUA=CreateMUA(CtrlVar,connectivity,coordinates);

if  CtrlVar.doplots==1 && CtrlVar.PlotMesh==1
     figure(1500+CtrlVar.GmshMeshingAlgorithm) ; hold off ; PlotFEmesh(coordinates,connectivity,CtrlVar)
end




end

