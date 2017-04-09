
function [UserVar,RunInfo,MUAnew,BCsNew,Fnew,lnew,GFnew]=AdaptMesh(UserVar,RunInfo,CtrlVar,MUAold,BCsOld,Fold,lold,GFold,RuvOld,Lubvb)

%save TestSave ; error('sdfa')

narginchk(10,10)
nargoutchk(7,7)

persistent FigMesh VideoMesh


MUAnew=MUAold;
Fnew=Fold;
BCsNew=BCsOld;
lnew=lold;
GFnew=GFold;


if CtrlVar.InfoLevelAdaptiveMeshing>=1
    fprintf('Before remeshing: '); PrintInfoAboutElementsSizes(CtrlVar,MUAold)
end


% I want to control directly when FE-mesh is plotted, and avoid the potential FE-mesh plotting in genmesh2d




if CtrlVar.doplots && CtrlVar.doAdaptMeshPlots && CtrlVar.InfoLevelAdaptiveMeshing>=2
    figure
    PlotFEmesh(MUAold.coordinates,MUAold.connectivity,CtrlVar);
    title(sprintf('Initial mesh on input to AdaptMesh : #Ele=%-i, #Nodes=%-i, #nod=%-i',MUAold.Nele,MUAold.Nnodes,MUAold.nod))
end


% Either there is mesh advance and/or retreat, or mesh adaptation
% If MeshAdvanceRetreat, do not do MeshAdapt

isMeshAdvanceRetreat = CtrlVar.FEmeshAdvanceRetreat && ( ReminderFraction(CtrlVar.time,CtrlVar.FEmeshAdvanceRetreatDT)<1e-5 || CtrlVar.FEmeshAdvanceRetreatDT==0);

isMeshAdapt=CtrlVar.AdaptMesh  ...
    && (CtrlVar.AdaptMeshInitial || (mod(CtrlVar.CurrentRunStepNumber,CtrlVar.AdaptMeshInterval)==0 ...
    && CtrlVar.AdaptMeshInterval>0 && CtrlVar.CurrentRunStepNumber>1)) ...
    && ~isMeshAdvanceRetreat;



if isMeshAdvanceRetreat
    
    [UserVar,RunInfo,MUAnew]=MeshAdvanceRetreat(UserVar,RunInfo,CtrlVar,MUAold,BCsOld,Fold,lold,GFold,RuvOld,Lubvb);
    
    if MUAnew.Nele==0
        fprintf('No elements left in mesh! \n ')
        return
    end
    
elseif isMeshAdapt
    
    for JJ=1:CtrlVar.AdaptMeshIterations
        
        MUAold=MUAnew;
        Fold=Fnew;
        BCsOld=BCsNew;
        GFold=GFnew;
        
        if CtrlVar.InfoLevelAdaptiveMeshing>=1
            fprintf(CtrlVar.fidlog,' =====  Remeshing at start of run step %-i. Iteration #%-i out of %-i \n ',CtrlVar.CurrentRunStepNumber,JJ,CtrlVar.AdaptMeshIterations);
        end
        
        %
        %  1) Determine new desired element sizes and identify elements for refinement
        %  or coarsening. 
        %
        %  2) Remesh: either global-remeshing or local-mesh refinement.
        %
        %  3) Map all variables onto new mesh.
        %
        
        
        [UserVar,RunInfo,xNod,yNod,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened]=...
            NewDesiredEleSizesAndElementsToRefineOrCoarsen(UserVar,RunInfo,CtrlVar,MUAold,BCsOld,Fold,lold,GFold,RuvOld,Lubvb);
        
        [UserVar,RunInfo,CtrlVar,MUAnew]=...
            Remeshing(UserVar,RunInfo,CtrlVar,MUAold,BCsOld,Fold,lold,GFold,...
            xNod,yNod,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened);
        
        
        
        if MUAnew.Nele==0
            fprintf('No elements left in mesh! \n ')
            return
        end
        
        
        [UserVar,Fnew,BCsNew,GFnew,lnew]=MapFbetweenMeshes(UserVar,CtrlVar,MUAold,MUAnew,Fold,BCsOld,GFold,lold);
        
        
        if (CtrlVar.doplots && CtrlVar.doAdaptMeshPlots && CtrlVar.InfoLevelAdaptiveMeshing>=1) || CtrlVar.CreateMeshAdaptVideo
            
            
            if isempty(FigMesh)
                FigMesh=figure ;
                %FigMesh.Position=[0 0 figsWidth 3*figHeights];
                FigMesh.Position=[1 1 1000 1000];% full laptop window
                axis tight
                if CtrlVar.CreateMeshAdaptVideo
                    VideoMesh=VideoWriter('MeshAdapt.avi');
                    VideoMesh.FrameRate=1;
                    open(VideoMesh);
                end
            else
                FigMesh=figure ;
                hold off
            end

            
            
            if CtrlVar.PlotBCs
                PlotBoundaryConditions(CtrlVar,MUAnew,BCsNew);
            else
                
                PlotFEmesh(MUAnew.coordinates,MUAnew.connectivity,CtrlVar); hold on
                title(sprintf('AdaptMesh Iteration %i: #Ele=%-i, #Nodes=%-i, #nod=%-i',JJ,MUAnew.Nele,MUAnew.Nnodes,MUAnew.nod))
            end
            
            if CtrlVar.GLmeshing
                hold on ; plot(xGLmesh,yGLmesh,'-rx','LineWidth',2);
            end
            
            axis tight
            if CtrlVar.CreateMeshAdaptVideo
                frame = getframe(gcf);
                writeVideo(VideoMesh,frame);
            end
        end
        
    end
    
end



if CtrlVar.InfoLevelAdaptiveMeshing>=1
    fprintf('After remeshing: ') ; 
    PrintInfoAboutElementsSizes(CtrlVar,MUAnew)
end


if CtrlVar.AdaptMeshAndThenStop
    
    if CtrlVar.doplots  && CtrlVar.InfoLevelAdaptiveMeshing>=10
        
        %%[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,varargin)
        xGL=[] ; yGL=[]; GLgeo=[];
        CtrlVar.PlotGLs=1;
        figure ; PlotMuaMesh(CtrlVar,MUAnew,[],CtrlVar.MeshColor);
        title(sprintf('Mesh after remeshing  \t #Ele=%-i, #Nodes=%-i, #nod=%-i',MUAnew.Nele,MUAnew.Nnodes,MUAnew.nod))
        hold on ;  [xGL,yGL]=PlotGroundingLines(CtrlVar,MUAnew,GFnew,GLgeo,xGL,yGL,'r');
        
        xGL=[] ; yGL=[]; GLgeo=[];
        figure ; PlotMuaMesh(CtrlVar,MUAold,[],CtrlVar.MeshColor);
        hold on ;  [xGL,yGL]=PlotGroundingLines(CtrlVar,MUAold,GFold,GLgeo,xGL,yGL,'r');
        title(sprintf('Mesh before remeshing  \t #Ele=%-i, #Nodes=%-i, #nod=%-i',MUAold.Nele,MUAold.Nnodes,MUAold.nod))
        %%
    end
    
    MUA=MUAnew;
    save(CtrlVar.SaveInitialMeshFileName,'MUA') ;
    fprintf(CtrlVar.fidlog,'New mesh was saved in %s .\n',CtrlVar.SaveAdaptMeshFileName);
    fprintf('Exiting after remeshing because CtrlVar.AdaptMeshAndThenStop set to true. \n ')
    return
end



%% MUA change, if so map all variables over to new mesh


[UserVar,Fnew,BCsNew,GFnew,lnew]=MapFbetweenMeshes(UserVar,CtrlVar,MUAold,MUAnew,Fold,BCsOld,GFold,lold);

% I still need to determine when to recalcuate uv 
% [UserVar,RunInfo,Fnew,lnew]= uv(UserVar,RunInfo,CtrlVar,MUAnew,BCsNew,Fnew,lnew);  % should really not be needed



if ~isempty(CtrlVar.SaveAdaptMeshFileName)
    MUA=MUAnew;
    save(CtrlVar.SaveAdaptMeshFileName,'MUA')
    fprintf(' Adapted FE mesh was saved in %s .\n',CtrlVar.SaveAdaptMeshFileName);
end

if ~isempty(VideoMesh)
    close(VideoMesh)
end


end

