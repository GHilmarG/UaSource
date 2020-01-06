
function [UserVar,RunInfo,MUAnew,BCsNew,Fnew,lnew]=AdaptMesh(UserVar,RunInfo,CtrlVar,MUAold,BCsOld,Fold,lold,RuvOld,Lubvb)


persistent AdaptMeshTime

narginchk(9,9)
nargoutchk(6,6)

%% Do all mesh modifications on MUAnew
% Only use MUAold and Fold for mapping of variables to the new mesh, ie always
% to the mapping from the original input mesh using the original input field
% values.

MUAnew=MUAold;
Fnew=Fold;
BCsNew=BCsOld;
lnew=lold;
RuvNew=RuvOld;

%%
isNewOutsideNodes=false  ; % true if during remeshing, in particular during manual deactivation of eliments, 
                           % some new nodes are introduced that are OUTSIDE of any of the elements in the old mesh

%%
xGLold=[] ; yGLold=[]; GLgeoold=[];
%%


isMeshAdvanceRetreat = CtrlVar.FEmeshAdvanceRetreat && ( ReminderFraction(CtrlVar.time,CtrlVar.FEmeshAdvanceRetreatDT)<1e-5 || CtrlVar.FEmeshAdvanceRetreatDT==0);

% isMeshAdapt=CtrlVar.AdaptMesh  ...
%     && (CtrlVar.AdaptMeshInitial || (mod(CtrlVar.CurrentRunStepNumber,CtrlVar.AdaptMeshRunStepInterval)==0 ...
%     && CtrlVar.AdaptMeshRunStepInterval>0 && CtrlVar.CurrentRunStepNumber>1)) ...
%     && ~isMeshAdvanceRetreat;

%%


if isempty(AdaptMeshTime)
    if CtrlVar.AdaptMeshTimeInterval>0
        AdaptMeshTime=ceil(CtrlVar.time/CtrlVar.AdaptMeshTimeInterval)*CtrlVar.AdaptMeshTimeInterval;
    else
        AdaptMeshTime=0;
    end
end

if CtrlVar.AdaptMeshTimeInterval==0
    isAdaptMeshTime=true;
elseif CtrlVar.time >= AdaptMeshTime
    isAdaptMeshTime=true;
    AdaptMeshTime=ceil((CtrlVar.time+eps)/CtrlVar.AdaptMeshTimeInterval)*CtrlVar.AdaptMeshTimeInterval;
else
    isAdaptMeshTime=false;
end

if CtrlVar.AdaptMeshRunStepInterval==0
    isAdaptMeshRunStepInterval=true;
elseif mod(CtrlVar.CurrentRunStepNumber,CtrlVar.AdaptMeshRunStepInterval)==0
    isAdaptMeshRunStepInterval=true;
else
    isAdaptMeshRunStepInterval=false;
end


%%

isMeshAdapt=CtrlVar.AdaptMesh  ...
    && (CtrlVar.AdaptMeshInitial || (isAdaptMeshRunStepInterval && isAdaptMeshTime)) ...
    && ~isMeshAdvanceRetreat;



if ~isMeshAdapt && ~isMeshAdvanceRetreat && ~CtrlVar.ManuallyDeactivateElements
    return
end



JJ=0 ;
nNewElements=inf;
nNewNodes=inf;
nNoChange=0;
RunInfo.MeshAdapt.isChanged=false;


CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;

if CtrlVar.InfoLevelAdaptiveMeshing>=1
    fprintf('Before remeshing: '); PrintInfoAboutElementsSizes(CtrlVar,MUAold);
end

OutsideValue=DefineOutsideValues(UserVar,CtrlVar,MUAold,Fold);


if isMeshAdvanceRetreat ||  isMeshAdapt
    
    if isMeshAdvanceRetreat
        
        [UserVar,RunInfo,MUAnew]=MeshAdvanceRetreat(UserVar,RunInfo,CtrlVar,MUAold,BCsOld,Fold,lold,Fold.GF,RuvNew,Lubvb);
        
        if MUAnew.Nele==0
            fprintf('No elements left in mesh! \n ')
            return
        end
        
    elseif isMeshAdapt
        
        while true
            
            JJ=JJ+1;
            
            if JJ>CtrlVar.AdaptMeshMaxIterations
                if CtrlVar.InfoLevelAdaptiveMeshing>=1
                    fprintf('Breaking out of adapt mesh iteration because number of iterations (%i) exceeds ''CtrlVar.AdaptMeshMaxIterations'' (%i)\n',...
                        JJ,CtrlVar.AdaptMeshMaxIterations)
                end
                break
            end
            
            if abs(nNewElements)<=CtrlVar.AdaptMeshUntilChangeInNumberOfElementsLessThan
                
                nNoChange=nNoChange+1;
                %if ~contains(RunInfo.MeshAdapt,'Bisection Coarsening')
                if nNoChange>1 || (nNewElements==0 && nNewNodes==0)
                    if CtrlVar.InfoLevelAdaptiveMeshing>=1
                        fprintf('Breaking out of adapt mesh iteration because change in the number of elements (%i) less than ''CtrlVar.AdaptMeshUntilChangeInNumberOfElementsLessThan'' (%i)\n',...
                            abs(nNewElements),CtrlVar.AdaptMeshUntilChangeInNumberOfElementsLessThan)
                    end
                    break
                end
            else
                nNoChange=0;
            end
            
            if CtrlVar.InfoLevelAdaptiveMeshing>=1
                fprintf(CtrlVar.fidlog,' --------->  Remeshing at start of run step %-i. Remeshing iteration #%-i (#Ele=%i,#Nodes=%i) \n ',CtrlVar.CurrentRunStepNumber,JJ,MUAnew.Nele,MUAnew.Nnodes);
            end
            
            %  Determine new desired element sizes and identify elements for refinement
            %  or coarsening.
            %  Note: F will be different on return if a new uv calculation
            %  needs to be done 
            [UserVar,RunInfo,Fnew,xNod,yNod,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened]=...
                NewDesiredEleSizesAndElementsToRefineOrCoarsen2(UserVar,RunInfo,CtrlVar,MUAnew,BCsNew,Fnew,lnew,Fnew.GF,RuvNew,Lubvb);
            
            %  Remesh: either global-remeshing or local-mesh refinement.
            %
            % If a local refinement/unrefinement is done using the newest-vertec bisection method,
            % the original MUA.RefineMesh structure is used. If the number of elements
            % changes, MUA is recreated using the elements and the coordinates in
            % MUA.RefineMesh.  Therefore, if some elements within MUA where previously
            % deactivated, the 'mother' elements will be reintroduced. Hence, I need to make
            % sure that any such elements are not deactivated again if manual deactivation
            % option is being used
            %
            
            NeleBefore=MUAnew.Nele;
            NnodesBefore=MUAnew.Nnodes;
            [UserVar,RunInfo,CtrlVar,MUAnew]=...
                Remeshing(UserVar,RunInfo,CtrlVar,MUAnew,BCsNew,Fnew,lnew,Fnew.GF,...
                xNod,yNod,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened);
            
            % if MUA changed check here if elements need to be deactivated
            
            if MUAnew.Nele==0
                fprintf('No elements left in mesh! \n ')
                return
            end
            
            nNewElements=MUAnew.Nele-NeleBefore;
            nNewNodes=MUAnew.Nnodes-NnodesBefore;
                        
            [UserVar,RunInfo,Fnew,BCsNew,lnew]=MapFbetweenMeshes(UserVar,RunInfo,CtrlVar,MUAold,MUAnew,Fold,BCsOld,lold,OutsideValue);
            
            if RunInfo.Mapping.nNotIdenticalNodesOutside>0
                isNewOutsideNodes=true  ; % true if during remeshing, in particular during manual deactivation of eliments, 
            end
            
            
            
            %% Plots
            if  CtrlVar.doplots && CtrlVar.doAdaptMeshPlots && CtrlVar.InfoLevelAdaptiveMeshing>=100
                
                FigureName='Adapt Mesh: before and after'; 
                fig=FindOrCreateFigure(FigureName);
                clf(fig);
                subplot(2,1,1)
                hold off
                
                PlotMuaMesh(CtrlVar,MUAold,[],CtrlVar.MeshColor);
                hold on ;  [xGLold,yGLold]=PlotGroundingLines(CtrlVar,MUAold,Fold.GF,GLgeoold,xGLold,yGLold,'r');
                title(sprintf('Before remeshing  \t #Ele=%-i, #Nodes=%-i, #nod=%-i',MUAold.Nele,MUAold.Nnodes,MUAold.nod))
                axis tight
                
                subplot(2,1,2)
                hold off
                xGL=[] ; yGL=[]; GLgeo=[];
                CtrlVar.PlotGLs=1;
                PlotMuaMesh(CtrlVar,MUAnew,[],CtrlVar.MeshColor);
                title(sprintf('After remeshing iteration #%i \t #Ele=%-i, #Nodes=%-i, #nod=%-i \n Change in the numbers of ele and nodes in current iteration is %i and %i ',...
                    JJ,MUAnew.Nele,MUAnew.Nnodes,MUAnew.nod,nNewElements,nNewNodes))
                hold on ;  [xGL,yGL]=PlotGroundingLines(CtrlVar,MUAnew,Fnew.GF,GLgeo,xGL,yGL,'r');
                axis tight
                
                fig.Children(2).XLim=fig.Children(1).XLim;
                fig.Children(2).YLim=fig.Children(1).YLim;
                sgtitle(sprintf('Adapt meshing at runstep %-i and time %f',CtrlVar.CurrentRunStepNumber,CtrlVar.time))
                
            end
            
        end
    end
    
    MUAnew=UpdateMUA(CtrlVar,MUAnew);
    
    % I now have done either MeshAdvanceRetreat or MeshAdapt, and I now update the
    % old MUA to the new MUA ahead of any further possible mesh modifications.
    OutsideValue.h=2*CtrlVar.ThickMin ;
    OutsideValue.s=mean(Fold.S)+CtrlVar.ThickMin*(1-mean(Fold.rho)/Fold.rhow);
    OutsideValue.b=OutsideValue.s-OutsideValue.h;
    OutsideValue.ub=0;
    OutsideValue.vb=0;
    
    [UserVar,RunInfo,Fnew,BCsNew,lnew]=MapFbetweenMeshes(UserVar,RunInfo,CtrlVar,MUAold,MUAnew,Fold,BCsOld,lold,OutsideValue);

end

%%

if CtrlVar.ManuallyDeactivateElements
    
    
    if CtrlVar.InfoLevelAdaptiveMeshing>=1
        fprintf('Manual deactivation of elements.\n');
    end
    
    % I need to have saved the original mesh if I want to be able to reactivate
    % regions. If the newest-vertex-biscetion local mesh-refinement option is being
    % used, then I can use MUA.RefineMesh for this purpose. However, if no such
    % local mesh refinement is used, I need to create this structure here.
    
    if  ~isfield(MUAnew,'RefineMesh')  ||  isempty(MUAnew.RefineMesh)
        mesh = genMesh(MUAnew.connectivity, MUAnew.coordinates);
        mesh.bd=[];
        mesh = genBisectionMesh(mesh);
        mesh = SelectRefinementEdge(mesh);
        if MUAold.nod~=3
            mesh.TR=triangulation(mesh.elements,mesh.coordinates);
        end
        MUAnew.RefineMesh=mesh;
    end
    
    
    % If the user wants to manually deactivate elements, I must re-introduce the full mesh
    % at each mesh refinement stage to allow for the re-activation of regions previously
    % deactivated. As this is also done within LocalMeshRefinement using newest-vertex
    % bisection in combination with manual deactivation of elements, this is only
    % occationally required.
    if  ~(size(MUAnew.RefineMesh.elements,1)==MUAnew.Nele && size(MUAnew.RefineMesh.coordinates,1)==MUAnew.Nnodes)

        MUAnew=CreateMUA(CtrlVar,MUAnew.RefineMesh.elements,MUAnew.RefineMesh.coordinates,MUAnew.RefineMesh);
        % The user might need estimates over the full mesh when making decisions, hence a
        % mapping to the new (full domain) mesh ahead of a call to
        % DefineElementsToDeactivate.m
        % it's enough to do this here because the mapping is otherwise always done in the Remeshing
        [UserVar,RunInfo,Fnew,BCsNew,lnew]=MapFbetweenMeshes(UserVar,RunInfo,CtrlVar,MUAold,MUAnew,Fold,BCsOld,lold,OutsideValue);
    
    end

    
    
    ElementsToBeDeactivated=false(MUAnew.Nele,1);
    
    [UserVar,ElementsToBeDeactivated]=...
        DefineElementsToDeactivate(UserVar,RunInfo,CtrlVar,MUAnew,MUAnew.xEle,MUAnew.yEle,ElementsToBeDeactivated,Fnew.s,Fnew.b,Fnew.S,Fnew.B,Fnew.rho,Fnew.rhow,Fnew.ub,Fnew.vb,Fnew.ud,Fnew.vd,Fnew.GF);
    
    if CtrlVar.doplots && CtrlVar.doAdaptMeshPlots && CtrlVar.InfoLevelAdaptiveMeshing>=100
        
        FigureName="Elements to be deactivated (red)";  
        fig=FindOrCreateFigure(FigureName);
        
        PlotMuaMesh(CtrlVar,MUAnew,ElementsToBeDeactivated,'r');
        hold on
        PlotMuaMesh(CtrlVar,MUAnew,~ElementsToBeDeactivated,'k');
        title('Elements to be deactivated in red')
    end
    
    MUAnew=DeactivateMUAelements(CtrlVar,MUAnew,ElementsToBeDeactivated);
    MUAnew=UpdateMUA(CtrlVar,MUAnew);
    
    CtrlVar.MapOldToNew.method="scatteredInterpolant"; % testing
    %
    
    OutsideValue.h=CtrlVar.ThickMin ;
    OutsideValue.s=mean(Fold.S)+CtrlVar.ThickMin*(1-mean(Fold.rho)/Fold.rhow);
    OutsideValue.b=OutsideValue.s-OutsideValue.h;
    OutsideValue.ub=0;
    OutsideValue.vb=0;
    [UserVar,RunInfo,Fnew,BCsNew,lnew]=MapFbetweenMeshes(UserVar,RunInfo,CtrlVar,MUAold,MUAnew,Fold,BCsOld,lold,OutsideValue);
    
    if RunInfo.Mapping.nNotIdenticalNodesOutside>0
        isNewOutsideNodes=true  ; % true if during remeshing, in particular during manual deactivation of eliments,
    end
end
%%


if CtrlVar.InfoLevelAdaptiveMeshing>=1
    fprintf('After remeshing: ') ;
    PrintInfoAboutElementsSizes(CtrlVar,MUAnew) ;
end

%% Before and after plot
%% Plots
if  CtrlVar.doplots && CtrlVar.doAdaptMeshPlots && CtrlVar.InfoLevelAdaptiveMeshing>=10
    
    FigureName='Adapt Mesh: before and after'; 
    fig=FindOrCreateFigure(FigureName);
    
    subplot(2,1,1)
    hold off
    xGL=[] ; yGL=[]; GLgeo=[];
    PlotMuaMesh(CtrlVar,MUAold,[],CtrlVar.MeshColor);
    hold on ;  [xGL,yGL]=PlotGroundingLines(CtrlVar,MUAold,Fold.GF,GLgeo,xGL,yGL,'r');
    title(sprintf('Before remeshing \t #Ele=%-i, #Nodes=%-i, #nod=%-i',MUAold.Nele,MUAold.Nnodes,MUAold.nod))
    axis tight
    
    subplot(2,1,2)
    hold off
    xGL=[] ; yGL=[]; GLgeo=[];
    CtrlVar.PlotGLs=1;
    PlotMuaMesh(CtrlVar,MUAnew,[],CtrlVar.MeshColor);
    title(sprintf('After remeshing  \t #Ele=%-i, #Nodes=%-i, #nod=%-i',MUAnew.Nele,MUAnew.Nnodes,MUAnew.nod))
    hold on ;  [xGL,yGL]=PlotGroundingLines(CtrlVar,MUAnew,Fnew.GF,GLgeo,xGL,yGL,'r');
    axis tight

     sgtitle(sprintf('Adapt meshing at runstep %-i and time %f',CtrlVar.CurrentRunStepNumber,CtrlVar.time))
end


if ~isempty(CtrlVar.SaveAdaptMeshFileName)
    MUA=MUAnew;
    save(CtrlVar.SaveAdaptMeshFileName,'MUA') ;
    fprintf(CtrlVar.fidlog,'New mesh was saved in %s .\n',CtrlVar.SaveAdaptMeshFileName);
end

if CtrlVar.AdaptMeshAndThenStop
    return
end



%%
%  Do velocities need to be recalculated?
%
%  Always recalculate velocities if:
%   CtrlVar.InitialDiagnosticStepAfterRemeshing is true
%   but also if mesh refinement method was not 'newest vertex bisection'
%
isMeshingLocalWithoutSmoothing=(contains(CtrlVar.MeshRefinementMethod,"explicit:local:red-green","IgnoreCase",true) && CtrlVar.LocalAdaptMeshSmoothingIterations==0) ...
    || contains(CtrlVar.MeshRefinementMethod,"local:newest vertex bisection","IgnoreCase",true);

isMeshChanged=HasMeshChanged(MUAold,MUAnew);

isRecalculateVelocities=isMeshChanged ||  isNewOutsideNodes  || CtrlVar.InitialDiagnosticStepAfterRemeshing || ~isMeshingLocalWithoutSmoothing ;

if ~CtrlVar.AdaptMeshAndThenStop
    if isRecalculateVelocities
        [UserVar,RunInfo,Fnew,lnew]= uv(UserVar,RunInfo,CtrlVar,MUAnew,BCsNew,Fnew,lnew);
    end
end

if ~isfield(RunInfo.MeshAdapt.Mesh,'Nele')
    RunInfo.MeshAdapt.Mesh.Nele=NaN;
    RunInfo.MeshAdapt.Mesh.Nnodes=NaN;
    RunInfo.MeshAdapt.Mesh.RunStepNumber=NaN;
    RunInfo.MeshAdapt.Mesh.time=NaN;
end

k=find(isnan(RunInfo.MeshAdapt.Mesh.Nele),1);

if isempty(k)
    RunInfo.MeshAdapt.Mesh.Nele=[RunInfo.MeshAdapt.Mesh.Nele;RunInfo.MeshAdapt.Mesh.Nele*0+NaN];
    RunInfo.MeshAdapt.Mesh.Nnodes=[RunInfo.MeshAdapt.Mesh.Nnodes;RunInfo.MeshAdapt.Mesh.Nnodes*0+NaN];
    RunInfo.MeshAdapt.Mesh.RunStepNumber=[RunInfo.MeshAdapt.Mesh.RunStepNumber;RunInfo.MeshAdapt.Mesh.RunStepNumber*0+NaN];
    RunInfo.MeshAdapt.Mesh.time=[RunInfo.MeshAdapt.Mesh.time;RunInfo.MeshAdapt.Mesh.time*0+NaN];
    
    k=find(isnan(RunInfo.MeshAdapt.Mesh.Nele),1);
end

RunInfo.MeshAdapt.Mesh.Nele(k)=MUAnew.Nele;
RunInfo.MeshAdapt.Mesh.Nnodes(k)=MUAnew.Nnodes;
RunInfo.MeshAdapt.Mesh.RunStepNumber(k)=CtrlVar.CurrentRunStepNumber;
RunInfo.MeshAdapt.Mesh.time(k)=CtrlVar.time;




%%

%%
%save TestSave
%[UserVar,RunInfo,F,l,BCs,GF,dt]=uvh(UserVar,RunInfo,CtrlVar,MUAnew,Fnew,Fnew,lnew,lnew,BCsNew);
%%



end

