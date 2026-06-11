function io = CreateGmshInitialMeshingInputFile(CtrlVar,MeshBoundaryCoordinates,FileName)

%
% io = CreateGmshInitialMeshingInputFile(CtrlVar,MeshBoundaryCoordinates,FileName)
% writes a .geo file which can be used as an input file for gmsh
%
% takes MeshBoundaryCoordinates which is a Npoints x 2 matrix with xy coordinates of the
% boundary of the domain.
% Also uses:
% CtrlVar.MeshSize  :  to define traget mesh size at boundary points
%
%
% The syntax for the MeshBoundaryCoordinates
%
%  [MeshId NaN ; x1 y1 ; x2 y2 ;  ... ; xn yn ; MeshId NaN ; x1 y1 ; x2 y2 ;  ... ; xn yn ;
%  MeshId identies to which mesh the line defined by the following x, y coordinates belongs
%  Individual lines are separted by the NaN following the MeshId
%  If there is only one line, there is not need to specify the MeshId
%  So if for example there is only one mesh surrounded by one line with no internal holes the syntax is simply
% [x1 y2 ; x2 y2 ; ... ; xn yn]
%
% The above syntax assumes that each line only belongs to one mesh.  If the same line belongs to
% more than one meshes then we need to specify
% the 'plane surface' (gmsh terminology) separatly. This is done in CtrlVar.GmshPlaneSurface
%
% Examples can be found in `ExamplesOFMeshGeneration.m'


io=0;  %  io=0 indicates no errors

FileName=[FileName,'.geo'] ;
fprintf(CtrlVar.fidlog,' Creating an input file for gmsh : %s \n',FileName) ;

[fileID,errmsh] = fopen(FileName,'w');

if fileID<0
   fprintf('opening the file %s resulted in an error:\n',FileName)
   disp(errmsh)
   error('Error opening a file. Possibly problems with permissions.')
end
    


if CtrlVar.GmshInputFormat==1
    
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
    
    % if no surface ids are given, set the surface id to 1 for all lines
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
        
        switch lower(CtrlVar.GmshBoundaryType)
            
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
                fprintf(CtrlVar.fidlog,'CtrlVar.GmshBoundaryType must be: [lines,spline] but is %s \n',CtrlVar.GmshBoundaryType);
                error('case fell through')
        end
        
    end
    
    
    if isfield(CtrlVar,'GmshPlaneSurface') && ~isempty(CtrlVar.GmshPlaneSurface)
        for I=1:numel(CtrlVar.GmshPlaneSurface)
            if ~isempty(CtrlVar.GmshPlaneSurface{I})
                str=strjoin(arrayfun(@(x) {sprintf('%i',x)},CtrlVar.GmshPlaneSurface{I}),',');
                fprintf(fileID,'Plane Surface(%i) = {%s};\n',I,str);
            end
        end
    else
        
        PlaneSurfaceIds=unique(iSurface);
        
        for PlaneSurfaceLabel=PlaneSurfaceIds(:,1)'
            PrintField=find(iSurface==PlaneSurfaceLabel)';
            str=strjoin(arrayfun(@(x) {sprintf('%i',x)},PrintField),',');
            fprintf(fileID,'Plane Surface(%i) = {%s};\n',PlaneSurfaceLabel,str);
        end
    end
    
elseif CtrlVar.GmshInputFormat==2

    
    % Points
    ItemsPerLine=20;
    Points=CtrlVar.Gmsh.Points;
    nPoints=size(Points,1);
    lPoints=1:nPoints;
    PrintField=[lPoints(:) Points(:,1) Points(:,2) zeros(nPoints,1) zeros(nPoints,1)+CtrlVar.MeshSize];
    fprintf(fileID,'Point(%i) = {%f,%f,%f,%f};\n',PrintField');
    
    % lines
    iLines=0 ; 
    for I=1:numel(CtrlVar.Gmsh.Lines)
        if isempty(CtrlVar.Gmsh.Lines{I})
            continue
        end
        
        Lines=CtrlVar.Gmsh.Lines{I};
        MultiLine{I}=[];
        for J=1:numel(Lines)-1
            fprintf(fileID,'Line(%i)={%i,%i};\n',iLines+1,Lines(J),Lines(J+1));
            MultiLine{I}=[MultiLine{I} ; iLines+1];
            iLines=iLines+1;
        end
    end
    
    % Loops
    for I=1:numel(CtrlVar.Gmsh.Loops)
        if isempty(CtrlVar.Gmsh.Loops{I})
            continue
        end
        
        Loop=CtrlVar.Gmsh.Loops{I};
        
        l=[];
        for J=1:numel(Loop)
            l= [l ; sign(Loop(J))*MultiLine{abs(Loop(J))}];
        end
        
        fprintf(fileID,'Line Loop(%i)={%i',I,l(1));
        for nn=2:numel(l)
            fprintf(fileID,',%i',l(nn));
            if mod(nn,ItemsPerLine)==0
                fprintf(fileID,'\n             ');
            end
        end
        fprintf(fileID,'};\n');
        
        
    end
    
    % Surfaces CtrlVar.Gmsh.PlaneSurfaces(1) = {1,2};
    for I=1:numel(CtrlVar.Gmsh.PlaneSurfaces)
        if isempty(CtrlVar.Gmsh.PlaneSurfaces{I})
            continue
        end
        
        Surface=CtrlVar.Gmsh.PlaneSurfaces{I};
        
        fprintf(fileID,'Plane Surface(%i)={%i',I,Surface(1));
        for nn=2:numel(Surface)
            fprintf(fileID,',%i',Surface(nn));
            if mod(nn,ItemsPerLine)==0
                fprintf(fileID,'\n                  ');
            end
        end
        fprintf(fileID,'};\n');
        
    end
    
end



fprintf(fileID,'Mesh.Algorithm = %i ; \n', CtrlVar.GmshMeshingAlgorithm);
fprintf(fileID,'Mesh.CharacteristicLengthMin = %f ; \n',CtrlVar.MeshSizeMin);
fprintf(fileID,'Mesh.CharacteristicLengthMax = %f ; \n',CtrlVar.MeshSizeMax);

if CtrlVar.GmshCharacteristicLengthExtendFromBoundary
    fprintf(fileID,'Mesh.CharacteristicLengthExtendFromBoundary = %i ; \n',1);
else
    fprintf(fileID,'Mesh.CharacteristicLengthExtendFromBoundary = %i ; \n',0);
end

if CtrlVar.GmshCharacteristicLengthFromCurvature
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

