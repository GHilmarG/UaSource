function io = CreateGmshInitialMeshingInputFile(CtrlVar,MeshBoundaryCoordinates,FileName)

%
% io = CreateGmshInitialMeshingInputFile(CtrlVar,MeshBoundaryCoordinates,FileName)
% writes a .geo file which can be used as an input file for gmesh
%
% takes MeshBoundaryCoordinates which is a Npoints x 2 matrix with xy coordinates of the
% boundary of the domain.
% Also uses:
% CtrlVar.MeshSize  :  to define traget mesh size at boundary points
%
%

io=0;  %  io=0 indicates no errors

FileName=[FileName,'.geo'] ;
fprintf(CtrlVar.fidlog,' Creating an input file for gmesh : %s \n',FileName) ;

fileID = fopen(FileName,'w');

% loops are separated by NaN. 
% What follows is simplified if I just add NaN at the beginning and end of MeshBoundaryCoordinates

if ~isnan(MeshBoundaryCoordinates(1,2))
    MeshBoundaryCoordinates=[1 NaN ; MeshBoundaryCoordinates];
end

if ~isnan(MeshBoundaryCoordinates(end,2))
    MeshBoundaryCoordinates=[MeshBoundaryCoordinates ; NaN NaN];
end
% how many loops do we have?
Ind=find(isnan(MeshBoundaryCoordinates(:,2))) ;

iSurface=MeshBoundaryCoordinates(Ind,1) ;
iSurface=iSurface(1:end-1);

if numel(iSurface)>1 && all(isnan(iSurface(2:end)))
    iSurface(:)=1;
end


nLoops=numel(Ind)-1;  % this is the number of loops
labelsEnd=0;
for I=1:nLoops

    Loop=zeros(Ind(I+1)-Ind(I)-1,2);
    Loop(:,1)=MeshBoundaryCoordinates(Ind(I)+1:Ind(I+1)-1,1);
    Loop(:,2)=MeshBoundaryCoordinates(Ind(I)+1:Ind(I+1)-1,2);
    
    nPoints=length(Loop);
    
    labelsStart=labelsEnd+1;
    labelsEnd=labelsStart+nPoints-1;
    labels=(labelsStart:labelsEnd)';
    PrintField=[labels Loop zeros(nPoints,1) zeros(nPoints,1)+CtrlVar.MeshSize]';
    
    fprintf(fileID,'Point(%i) = {%f,%f,%f,%f};\n',PrintField);
    
    switch lower(CtrlVar.GmeshBoundaryType)
        
        case 'lines'
            % lines
            
            PrintField=[labels labels circshift(labels,-1)]';
            fprintf(fileID,'Line(%i) = {%i,%i};\n',PrintField);
            
            % line loop
            fprintf(fileID,'Line Loop(%i) = {%i:%i};\n',I,labels(1),labels(end));
            
        case 'spline'
            fprintf(fileID,'Spline(%i) = {%i:%i,%i};\n',labels(end)+1,labels(1),labels(end),1);
            fprintf(fileID,'Line Loop(%i) = {%i};\n',I,labels(end)+1);
        otherwise
            fprintf(CtrlVar.fidlog,'CtrlVar.GmeshBoundaryType must be: [lines,spline] but is %s \n',CtrlVar.GmeshBoundaryType);
            error('case fell through')
    end
    
end


PlaneSurfaceIds=unique(iSurface);

for PlaneSurfaceLabel=PlaneSurfaceIds(:,1)'
    PrintField=find(iSurface==PlaneSurfaceLabel)';
    str=strjoin(arrayfun(@(x) {sprintf('%i',x)},PrintField),',');
    fprintf(fileID,'Plane Surface(%i) = {%s};\n',PlaneSurfaceLabel,str);
end

fprintf(fileID,'Mesh.Algorithm = %i ; \n', CtrlVar.GmeshMeshingAlgorithm);
fprintf(fileID,'Mesh.CharacteristicLengthMin = %f ; \n',CtrlVar.MeshSizeMin);
fprintf(fileID,'Mesh.CharacteristicLengthMax = %f ; \n',CtrlVar.MeshSizeMax);

if CtrlVar.GmeshCharacteristicLengthExtendFromBoundary
    fprintf(fileID,'Mesh.CharacteristicLengthExtendFromBoundary = %i ; \n',1);
else
    fprintf(fileID,'Mesh.CharacteristicLengthExtendFromBoundary = %i ; \n',0);
end

if CtrlVar.GmeshCharacteristicLengthFromCurvature
    fprintf(fileID,'Mesh.CharacteristicLengthFromCurvature = 1 ; \n');
else
    fprintf(fileID,'Mesh.CharacteristicLengthFromCurvature = 0 ; \n');
end

% now add any additional input lines to gmsh given by user
for I=1:length(CtrlVar.GmshGeoFileAdditionalInputLines)
    fprintf(fileID,' %s \n ',CtrlVar.GmshGeoFileAdditionalInputLines{I});
end

%fprintf(fileID,'Mesh.CharacteristicLengthFromPoints = 0 ;');

fclose(fileID);

end

