function [UserVar,RunInfo,MUA]=GlobalRemeshing(UserVar,RunInfo,CtrlVar,MUA,xNod,yNod,EleSizeDesired,F)
    
    
    
    if CtrlVar.InfoLevelAdaptiveMeshing>=1
        
        fprintf('\t Global adaptive remeshing at run-step %i and time=%f \n',CtrlVar.CurrentRunStepNumber,CtrlVar.time)
        
    end
    
    
    
    if any(isnan(EleSizeDesired)) ; save TestSave ; error('fdsa') ; end
    
    %if strcmp(CtrlVar.MeshGenerator,'gmsh')
    if norm([xNod-MUA.coordinates(:,1);yNod-MUA.coordinates(:,2)]) > 100*eps
        error('RemeshingBasedOnExplicitErrorEstimate:gmsh','When using gmsh desired ele sizes must be defined at nodes')
    end
    TRIxy0=TriFE(MUA.connectivity);  % this is the triangulation of the input FEmesh over which
    % the error estimation is performed
    %end
    
    
    
    
    
    %% Global remeshing
    
    
    
    EleSizeScalarField.xy=[xNod(:) yNod];
    EleSizeScalarField.EleSize=EleSizeDesired(:);
    EleSizeScalarField.TRI=TRIxy0;
    
    
    %     if CtrlVar.GLmeshing==1
    %         GF = GL2d(F.B,F.S,F.h,F.rhow,F.rho,F.MUA.connectivity,CtrlVar);  % Do I need this?, broken anyhow...
    %         [MeshBoundaryCooWithGLcoo,edge,face,xGLmesh,yGLmesh]=glLineEdgesFaces(GF,MUA.coordinates,MUA.connectivity,MeshBoundaryCoordinates,CtrlVar);
    %         [UserVar,MUA]=genmesh2d(UserVar,CtrlVar,MeshBoundaryCooWithGLcoo,edge,face);
    %     else
    
    [UserVar,MUA]=genmesh2d(UserVar,CtrlVar,F,EleSizeScalarField);
    
    %    end
    
    %figure ; PlotFEmesh(coordinates,connectivity,CtrlVar);
    
    
    
    % if the ratio of actual number of elements to desired number of elements
    % is either larger than `EleFactorU' or smaller than 'EleFactorDown' then MeshSizeMin is scaled
    % up or down by the factor 'EleFactor' and a new remeshing is done.
    % The maximum size of elements does not change and is always MeshSizeMax. Only MeshSizeMin is changed
    % together with the error chriteria, MeshSizeMax and the number of elements are the main variables,
    % affecting the mesh. The new value of MeshSizeMin is written out and in principle one should use that
    % value of MeshSizeMin as an input value in future runs of the same model
    
    %EleFactorUp=1.3 ; EleFactorDown=0.1;
    EleFactorUp=CtrlVar.MaxNumberOfElementsUpperLimitFactor;
    EleFactorDown=CtrlVar.MaxNumberOfElementsLowerLimitFactor;
    
    It=0;
    
    while (MUA.Nele>EleFactorUp*CtrlVar.MaxNumberOfElements ||  MUA.Nele<EleFactorDown*CtrlVar.MaxNumberOfElements ) && It<4
        
        ScalingFactor=sqrt(MUA.Nele/CtrlVar.MaxNumberOfElements);
        
        %EleSize=EleSize*ScalingFactor;
        % make sure that no element is greater than max element size
        %CtrlVar.MeshSize(EleSize>CtrlVar.MeshSizeMax,3)=CtrlVar.MeshSizeMax;
        
        % rescale elements sizes so that maximum size is still MeshSizeMax but min ele size is scaled either up or down
        maxE=max(EleSizeDesired);
        minE=min(EleSizeDesired);
        
        NewMinE=min([minE*ScalingFactor,0.9*CtrlVar.MeshSizeMax]);
        if maxE~=minE
            EleSizeDesired=NewMinE+(CtrlVar.MeshSizeMax-NewMinE)*(EleSizeDesired-minE)/(maxE-minE);
        else
            EleSizeDesired=NewMinE+EleSizeDesired*0;
        end
        
        It=It+1;
        
        if NewMinE ~=CtrlVar.MeshSizeMin
            warning('RmeshingBasedOnExplicitErrorEstimate:NewMeshSizeMin','CtrlVar.MeshSizeMin set to %f',NewMinE)
        end
        
        CtrlVar.MeshSizeMin=NewMinE;
        
        
        %    switch lower(CtrlVar.MeshGenerator)
        %        case 'mesh2d'
        %            CtrlVar.MeshSize=zeros(length(xNod),3);
        %            CtrlVar.MeshSize(:,1)=xNod ; CtrlVar.MeshSize(:,2)=yNod; CtrlVar.MeshSize(:,3)=EleSizeDesired;
        %        case 'gmsh'
        EleSizeScalarField.xy=[xNod(:) yNod(:)] ;
        EleSizeScalarField.EleSize=EleSizeDesired(:) ;
        EleSizeScalarField.TRI=TRIxy0 ;
        
        %            if any(isnan(EleSizeDesired)) ; error('fdsa') ; end
        %        otherwise
        %            error('Mesh generator not correctly defined. Define variable CtrlVar.MeshGenerator {mesh2d|gmsh} ')
        %    end
        
        if CtrlVar.InfoLevelAdaptiveMeshing>=1
            if MUA.Nele>EleFactorUp*CtrlVar.MaxNumberOfElements
                fprintf(CtrlVar.fidlog,'Nele=%-i > %-i , hence desired meshsize is scaled up. New MeshSizeMin is %-g \n',MUA.Nele,CtrlVar.MaxNumberOfElements,CtrlVar.MeshSizeMin);
            else
                fprintf(CtrlVar.fidlog,'Nele=%-i < %-i , hence desired meshsize is scaled down. New MeshSizeMin is %-g \n',MUA.Nele,CtrlVar.MaxNumberOfElements,CtrlVar.MeshSizeMin);
            end
        end
        
        %         if CtrlVar.GLmeshing==1
        %             [UserVar,MUA]=genmesh2d(UserVar,CtrlVar,MeshBoundaryCooWithGLcoo,edge,face);
        %         else
        [UserVar,MUA]=genmesh2d(UserVar,CtrlVar,F,EleSizeScalarField);
        %end
        
        if CtrlVar.InfoLevelAdaptiveMeshing>=1
            fprintf(CtrlVar.fidlog,'new Nele after rescaling is %-i and CtrlVar.MaxNumberOfElements is %-i  \n',MUA.Nele,CtrlVar.MaxNumberOfElements);
        end
    end
    
    
    
end