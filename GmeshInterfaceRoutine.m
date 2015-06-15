function [coordinates,connectivity]=GmeshInterfaceRoutine(CtrlVar,MeshBoundaryCoordinates,GmeshBackgroundScalarField)

%
%  [coordinates,connectivity]=GmeshInterfaceRoutine(CtrlVar,MeshBoundaryCoordinates,GmeshBackgroundScalarField)
%
% requires first two input arguments, and uses the following CtrlVar fields:
%
% CtrlVar.GmeshMeshingMode   = [create new gmesh .geo input file | mesh domain | load .msh]
%       these can be combined as for example: create new gmesh .geo input file and mesh domain and load .msh file
%                                         or: mesh domain and load .msh file
%
% CtrlVar.GmeshFile  : a .geo input file for gmsh
%                           if GmeshMeshingMode='mesh domain and load .msh file'  then CtrlVar.GmeshFile must be an existion .geo file
%                           if GmeshMeshingMode='create new gmesh .geo input file and mesh domain and load .msh file'
%                           then CtrlVar.GmeshFile is first created and then used as an input file for gmesh
%


% get rid of eventual file extension
[pathstr,name,ext]=fileparts(CtrlVar.GmeshFile);
CtrlVar.GmeshFile=[pathstr,name];

% set the path to gmsh
str=computer;
if  strcmp(str,'GLNXA64')     % Unix
    gmshRunString='gmsh ';
else
    GmeshHomeDirectory=getenv('GmeshHomeDirectory');    % Windows
    gmshRunString=fullfile(GmeshHomeDirectory,'gmsh.exe ');
end



if nargin<3 || isempty(GmeshBackgroundScalarField)  % mesh domain (without using background mesh file giving desired ele sizes)
    
    
    % possibly create a new .geo input file for GMESH. If not then it is assumed that
    % such a file already exists
    if ~isempty(strfind(lower(CtrlVar.GmeshMeshingMode),'create new gmesh .geo input file'))
        CreateGmshInitialMeshingInputFile(CtrlVar,MeshBoundaryCoordinates,CtrlVar.GmeshFile);
    end
    % read .geo input file for GMESH
    [fid,Message]=fopen([CtrlVar.GmeshFile,'.geo']);
    if fid==-1
        fprintf(CtrlVar.fidlog,'%s \n',Message);
        fprintf(CtrlVar.fidlog,'Error opening file: %s \n',[CtrlVar.GmeshFile,'.geo']);
        
        if exist(fullfile(cd,[CtrlVar.GmeshFile,'.geo']), 'file')  == 0
            fprintf(CtrlVar.fidlog,' File: %s does not exist \n',[CtrlVar.GmeshFile,'.geo']);
        end
        error(' file input error ')
    end
    
    if ~isempty(strfind(lower(CtrlVar.GmeshMeshingMode),'mesh domain'))
        
        RunString=[gmshRunString,CtrlVar.GmeshFile,'.geo -2 -v 1'];
        
    end
    
else  % remesh with a given scalar background field defining desired ele sizes
    
    fprintf(CtrlVar.fidlog,'Remesh using a background scalar field \n');
    CreateGmshInitialMeshingInputFile(CtrlVar,MeshBoundaryCoordinates,CtrlVar.GmeshFile);
    
    FileName=[CtrlVar.GmeshFile,'.pos'];
    fprintf(CtrlVar.fidlog,'Creating a Gmesh scalar post file %s \n',FileName);
    io = CreateGmshBackgroundScalarMesh(GmeshBackgroundScalarField.xy,GmeshBackgroundScalarField.TRI,GmeshBackgroundScalarField.EleSize,FileName);
    
    RunString=[gmshRunString,CtrlVar.GmeshFile,'.geo -bgm ',FileName,' -2 -v 1'];
    
    
end

% call gmesh
if ~isempty(strfind(lower(CtrlVar.GmeshMeshingMode),'mesh domain'))
    fprintf(CtrlVar.fidlog,'Calling gmesh with %s as a geo input file, and creating %s output file \n',[CtrlVar.GmeshFile,'.geo'],[CtrlVar.GmeshFile,'.msh']);
    fprintf(CtrlVar.fidlog,'The gmesh call is:  %s  \n',RunString);
    status=system(RunString);
    
    if status~=0
        fprintf(CtrlVar.fidlog,'gmesh returns with an error running (%s) \n',RunString);
        fprintf(CtrlVar.fidlog,' A possible reason is that the environmetal variable GmeshHomeDirectory is not set correctly\n');
        error(' Fix this somehow!')
    end
end

%%
if ~isempty(strfind(lower(CtrlVar.GmeshMeshingMode),'load .msh'))
    FileName=[CtrlVar.GmeshFile,'.msh'];
    fprintf(CtrlVar.fidlog,'Loading %s \n',FileName);
    Gmesh=load_gmshGHG(FileName);
    
    
    connectivity=Gmesh.TRIANGLES(1:Gmesh.nbTriangles,1:3);
    coordinates=Gmesh.POS(1:Gmesh.nbNod,1:2);
    
    if numel(coordinates)==0 ; error(' no coordinates in meshfile %s/%s \n',pwd,FileName) ; end
    
    [coordinates,connectivity]=RemoveNodesNotInConnectivity(coordinates,connectivity);
    
end

end



