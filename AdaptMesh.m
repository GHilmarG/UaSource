
function [UserVar,RunInfo,MUAnew,BCsNew,Fnew,lnew,GFnew]=AdaptMesh(UserVar,RunInfo,CtrlVar,MUAold,BCsOld,Fold,lold,GFold,RuvOld,Lubvb)

%save TestSave ; error('sdfa')

narginchk(10,10)
nargoutchk(7,7)

MUAnew=MUAold;
Fnew=Fold;
BCsNew=BCsOld;
lnew=lold;
GFnew=GFold;

CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;

if CtrlVar.InfoLevelAdaptiveMeshing>=1
    fprintf('Before remeshing: '); PrintInfoAboutElementsSizes(CtrlVar,MUAold)
end


isMeshAdvanceRetreat = CtrlVar.FEmeshAdvanceRetreat && ( ReminderFraction(CtrlVar.time,CtrlVar.FEmeshAdvanceRetreatDT)<1e-5 || CtrlVar.FEmeshAdvanceRetreatDT==0);

isMeshAdapt=CtrlVar.AdaptMesh  ...
    && (CtrlVar.AdaptMeshInitial || (mod(CtrlVar.CurrentRunStepNumber,CtrlVar.AdaptMeshInterval)==0 ...
    && CtrlVar.AdaptMeshInterval>0 && CtrlVar.CurrentRunStepNumber>1)) ...
    && ~isMeshAdvanceRetreat;


JJ=0 ;
nNewElements=inf;
nNewNodes=inf;
nNoChange=0;
RunInfo.MeshAdapt.isChanged=false;

if isMeshAdvanceRetreat
    
    [UserVar,RunInfo,MUAnew]=MeshAdvanceRetreat(UserVar,RunInfo,CtrlVar,MUAold,BCsOld,Fold,lold,GFold,RuvOld,Lubvb);
    
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
        
        
        
        MUAold=MUAnew;
        Fold=Fnew;
        BCsOld=BCsNew;
        GFold=GFnew;
        
        if CtrlVar.InfoLevelAdaptiveMeshing>=1
            fprintf(CtrlVar.fidlog,' =====  Remeshing at start of run step %-i. Remeshing iteration #%-i \n ',CtrlVar.CurrentRunStepNumber,JJ);
        end
        
        %  Determine new desired element sizes and identify elements for refinement
        %  or coarsening.
        [UserVar,RunInfo,Fold,xNod,yNod,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened]=...
            NewDesiredEleSizesAndElementsToRefineOrCoarsen2(UserVar,RunInfo,CtrlVar,MUAold,BCsOld,Fold,lold,GFold,RuvOld,Lubvb);
        
        %  Remesh: either global-remeshing or local-mesh refinement.
        [UserVar,RunInfo,CtrlVar,MUAnew]=...
            Remeshing(UserVar,RunInfo,CtrlVar,MUAold,BCsOld,Fold,lold,GFold,...
            xNod,yNod,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened);
        
        
        if MUAnew.Nele==0
            fprintf('No elements left in mesh! \n ')
            return
        end
        
        
        nNewElements=MUAnew.Nele-MUAold.Nele;
        nNewNodes=MUAnew.Nnodes-MUAold.Nnodes;
        
        [UserVar,RunInfo,Fnew,BCsNew,GFnew,lnew]=MapFbetweenMeshes(UserVar,RunInfo,CtrlVar,MUAold,MUAnew,Fold,BCsOld,GFold,lold);
        
        
        %% Plots
        if  CtrlVar.doplots && CtrlVar.doAdaptMeshPlots && CtrlVar.InfoLevelAdaptiveMeshing>=10
            
            FigNameBA='Adapt Mesh: before and after';
            fig=findobj(0,'name',FigNameBA);
            if isempty(fig)
                fig=figure('name',FigNameBA);
                fig.Position=[100,100,1000,1000] ;
            else
                fig=figure(fig);
                hold off
            end
            
            subplot(2,1,1)
            hold off
            xGL=[] ; yGL=[]; GLgeo=[];
            PlotMuaMesh(CtrlVar,MUAold,[],CtrlVar.MeshColor);
            hold on ;  [xGL,yGL]=PlotGroundingLines(CtrlVar,MUAold,GFold,GLgeo,xGL,yGL,'r');
            title(sprintf('Adapt-mesh iteration #%i \n Before remeshing  \t #Ele=%-i, #Nodes=%-i, #nod=%-i',JJ,MUAold.Nele,MUAold.Nnodes,MUAold.nod))
            axis tight
            
            subplot(2,1,2)
            hold off
            xGL=[] ; yGL=[]; GLgeo=[];
            CtrlVar.PlotGLs=1;
            PlotMuaMesh(CtrlVar,MUAnew,[],CtrlVar.MeshColor);
            title(sprintf('After remeshing  \t #Ele=%-i, #Nodes=%-i, #nod=%-i',MUAnew.Nele,MUAnew.Nnodes,MUAnew.nod))
            hold on ;  [xGL,yGL]=PlotGroundingLines(CtrlVar,MUAnew,GFnew,GLgeo,xGL,yGL,'r');
            axis tight
            
            
            
        end
        %%
    end
end




if ~isempty(CtrlVar.SaveAdaptMeshFileName)
    MUA=MUAnew;
    save(CtrlVar.SaveInitialMeshFileName,'MUA') ;
    fprintf(CtrlVar.fidlog,'New mesh was saved in %s .\n',CtrlVar.SaveAdaptMeshFileName);
end

%%
%% map variables to new mesh

[UserVar,RunInfo,Fnew,BCsNew,GFnew,lnew]=MapFbetweenMeshes(UserVar,RunInfo,CtrlVar,MUAold,MUAnew,Fold,BCsOld,GFold,lold);

%%

if CtrlVar.ManuallyDeactivateElements
    
    
    if CtrlVar.InfoLevelAdaptiveMeshing>=1
        fprintf('Manual deactivation of elements.\n'); 
    end
    
    
    xEle=Nodes2EleMean(MUAnew.connectivity,MUAnew.coordinates(:,1));
    yEle=Nodes2EleMean(MUAnew.connectivity,MUAnew.coordinates(:,2));
    ElementsToBeDeactivated=false(MUAnew.Nele,1);
    
    [UserVar,ElementsToBeDeactivated]=...
        DefineElementsToDeactivate(UserVar,RunInfo,CtrlVar,MUAnew,xEle,yEle,ElementsToBeDeactivated,Fnew.s,Fnew.b,Fnew.S,Fnew.B,Fnew.rho,Fnew.rhow,Fnew.ub,Fnew.vb,Fnew.ud,Fnew.vd,GFnew);
    
    if CtrlVar.doplots && CtrlVar.doAdaptMeshPlots
        figure
        PlotMuaMesh(CtrlVar,MUA,ElementsToBeDeactivated,'r')
        hold on
        PlotMuaMesh(CtrlVar,MUA,~ElementsToBeDeactivated,'k')
        title('Elements to be deactivated in red')
    end
    
    [MUAnew.coordinates,MUAnew.connectivity]=DeactivateElements(CtrlVar,ElementsToBeDeactivated,MUAnew.coordinates,MUAnew.connectivity);
    
    MUAnew=UpdateMUA(CtrlVar,MUAnew);
    %MUAnew=CreateMUA(CtrlVar,connectivity,coordinates);
    
    [UserVar,RunInfo,Fnew,BCsNew,GFnew,lnew]=MapFbetweenMeshes(UserVar,RunInfo,CtrlVar,MUAold,MUAnew,Fold,BCsOld,GFold,lold);
end
%%


if CtrlVar.InfoLevelAdaptiveMeshing>=1
    fprintf('After remeshing: ') ;
    PrintInfoAboutElementsSizes(CtrlVar,MUAnew)
end



if CtrlVar.AdaptMeshAndThenStop
    return
end



%%
%  Do velocities need to be recalculated?
%
%  Always recalculate velocities if:
%
%   CtrlVar.InitialDiagnosticStepAfterRemeshing is true
%   but also if mesh refinement method was not 'newest vertex bisection'
%
isMeshingLocalWithoutSmoothing=contains(CtrlVar.MeshRefinementMethod,'local','IgnoreCase',true) && CtrlVar.LocalAdaptMeshSmoothingIterations==0;

if ~CtrlVar.AdaptMeshAndThenStop
    if CtrlVar.InitialDiagnosticStepAfterRemeshing || ~isMeshingLocalWithoutSmoothing
        isMeshChanged=HasMeshChanged(MUAold,MUAnew);
        if isMeshChanged
            [UserVar,RunInfo,Fnew,lnew]= uv(UserVar,RunInfo,CtrlVar,MUAnew,BCsNew,Fnew,lnew);
        end
    end
end


end

