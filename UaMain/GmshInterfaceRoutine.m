function [coordinates,connectivity]=GmshInterfaceRoutine(CtrlVar,MeshBoundaryCoordinates,EleSizeScalarField)

%
%  [coordinates,connectivity]=GmshInterfaceRoutine(CtrlVar,MeshBoundaryCoordinates,EleSizeScalarField)
%
% requires first two input arguments, and uses the following CtrlVar fields:
%
% CtrlVar.GmshMeshingMode   = [create new gmsh .geo input file | mesh domain | load .msh]
%       these can be combined as for example: create new gmsh .geo input file and mesh domain and load .msh file
%                                         or: mesh domain and load .msh file
%
% CtrlVar.GmshFile  : a .geo input file for gmsh
%                           if GmshMeshingMode='mesh domain and load .msh file'  then CtrlVar.GmshFile must be an existion .geo file
%                           if GmshMeshingMode='create new gmsh .geo input file and mesh domain and load .msh file'
%                           then CtrlVar.GmshFile is first created and then used as an input file for gmsh
%


% get rid of eventual file extension
[pathstr,name]=fileparts(CtrlVar.GmshFile);
CtrlVar.GmshFile=[pathstr,name];

% set the path to gmsh

str=computer;
GmshHomeDirectory=getenv('GmshHomeDirectory');  
    
if  strcmp(str,'GLNXA64') || strcmp(str,'MACI64')
    gmshRunString=fullfile(GmshHomeDirectory,'gmsh ');     % Unix
else
    if isempty(GmshHomeDirectory)
        warning('Ua:GmshInterfaceRoutine','The environmental variable GmshHomeDirectory is empty')
    end
    gmshRunString=fullfile(GmshHomeDirectory,'gmsh.exe ');  
end



if nargin<3 || isempty(EleSizeScalarField)  % mesh domain (without using background mesh file giving desired ele sizes)
    
    
    % possibly create a new .geo input file for GMESH. If not then it is assumed that
    % such a file already exists
    if contains(lower(CtrlVar.GmshMeshingMode),'create new gmsh .geo input file')
        CreateGmshInitialMeshingInputFile(CtrlVar,MeshBoundaryCoordinates,CtrlVar.GmshFile);
    end
    % read .geo input file for GMESH
    [fid,Message]=fopen([CtrlVar.GmshFile,'.geo']);
    if fid==-1
        fprintf(CtrlVar.fidlog,'%s \n',Message);
        fprintf(CtrlVar.fidlog,'Error opening file: %s \n',[CtrlVar.GmshFile,'.geo']);
        
        if exist(fullfile(cd,[CtrlVar.GmshFile,'.geo']), 'file')  == 0
            fprintf(CtrlVar.fidlog,' File: %s does not exist \n',[CtrlVar.GmshFile,'.geo']);
        end
        error(' file input error ')
    end
    
    if contains(lower(CtrlVar.GmshMeshingMode),'mesh domain')
        
        RunString=[gmshRunString,CtrlVar.GmshFile,'.geo -2 -format msh2 -v ',num2str(CtrlVar.GmshVerbosityLevel)];
        
    end
    


else  % remesh with a given scalar background field defining desired ele sizes
    
    fprintf(CtrlVar.fidlog,'Remesh using a background scalar field \n');
    CreateGmshInitialMeshingInputFile(CtrlVar,MeshBoundaryCoordinates,CtrlVar.GmshFile);
    
    FileName=[CtrlVar.GmshFile,'.pos'];
    fprintf(CtrlVar.fidlog,'Creating a Gmsh scalar post file %s \n',FileName);
    CreateGmshBackgroundScalarMesh(EleSizeScalarField.xy,EleSizeScalarField.TRI,EleSizeScalarField.EleSize,FileName);
    
    RunString=[gmshRunString,CtrlVar.GmshFile,'.geo -bgm ',FileName,' -format msh2 -2 -v ',num2str(CtrlVar.GmshVerbosityLevel)];
    
    
end


% call gmsh
if contains(lower(CtrlVar.GmshMeshingMode),'mesh domain')
    fprintf(CtrlVar.fidlog,'Calling gmsh with %s as a geo input file, and creating %s output file \n',[CtrlVar.GmshFile,'.geo'],[CtrlVar.GmshFile,'.msh']);
    fprintf(CtrlVar.fidlog,'The gmsh call is:  %s  \n',RunString);
    if CtrlVar.GmshPause>0
        fprintf('Introducing a pause of %f sec before running gmsh.\n',CtrlVar.GmshPause)
        pause(CtrlVar.GmshPause);
        
    end
    status=system(RunString);
    
    if status~=0
        fprintf(CtrlVar.fidlog,'gmsh returns with an error running (%s) \n',RunString);
        fprintf(CtrlVar.fidlog,' A possible reason is that the environment variable GmshHomeDirectory is not set correctly\n');
        error('Ua:GmshInterfaceRoutine:GmshRunFailed',' Fix this somehow.')
    end
end

%%
if contains(lower(CtrlVar.GmshMeshingMode),'load .msh')
    FileName=[CtrlVar.GmshFile,'.msh'];
    fprintf(CtrlVar.fidlog,'Loading %s \n',FileName);
    Gmsh=load_gmshGHG(FileName);
    
    
    connectivity=Gmsh.TRIANGLES(1:Gmsh.nbTriangles,1:3);
    coordinates=Gmsh.POS(1:Gmsh.nbNod,1:2);
    
    if numel(coordinates)==0 
        error('Ua:GmshInterfaceRoutine:NoCoordinatesInMeshFile',' no coordinates in meshfile %s/%s \n',pwd,FileName) ;
    end
    
    [coordinates,connectivity]=RemoveNodesNotInConnectivity(coordinates,connectivity);
    
end

end



