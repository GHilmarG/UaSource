
function [UserVar,RunInfo,MUAnew,BCsNew,Fnew,lnew]=AdaptMesh(UserVar,RunInfo,CtrlVar,MUAold,BCsOld,Fold,lold,RuvOld,Lubvb)


persistent AdaptMeshTime


narginchk(9,9)
nargoutchk(6,6)

MUAnew=MUAold;
Fnew=Fold;
BCsNew=BCsOld;
lnew=lold;



%%


isMeshAdvanceRetreat = CtrlVar.FEmeshAdvanceRetreat && ( ReminderFraction(CtrlVar.time,CtrlVar.FEmeshAdvanceRetreatDT)<1e-5 || CtrlVar.FEmeshAdvanceRetreatDT==0);

% isMeshAdapt=CtrlVar.AdaptMesh  ...
%     && (CtrlVar.AdaptMeshInitial || (mod(CtrlVar.CurrentRunStepNumber,CtrlVar.AdaptMeshRunStepInterval)==0 ...
%     && CtrlVar.AdaptMeshRunStepInterval>0 && CtrlVar.CurrentRunStepNumber>1)) ...
%     && ~isMeshAdvanceRetreat;


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


isMeshAdapt=CtrlVar.AdaptMesh  ...
    && (CtrlVar.AdaptMeshInitial || (isAdaptMeshRunStepInterval && isAdaptMeshTime)) ...
    && ~isMeshAdvanceRetreat;


if ~isMeshAdapt && ~isMeshAdvanceRetreat
    return
end


JJ=0 ;
nNewElements=inf;
nNewNodes=inf;
nNoChange=0;
RunInfo.MeshAdapt.isChanged=false;


CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;

if CtrlVar.InfoLevelAdaptiveMeshing>=1
    fprintf('Before remeshing: '); PrintInfoAboutElementsSizes(CtrlVar,MUAold)
end



if isMeshAdvanceRetreat
    
    [UserVar,RunInfo,MUAnew]=MeshAdvanceRetreat(UserVar,RunInfo,CtrlVar,MUAold,BCsOld,Fold,lold,Fold.GF,RuvOld,Lubvb);
    
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
        
      
        
        if CtrlVar.InfoLevelAdaptiveMeshing>=1
            fprintf(CtrlVar.fidlog,' --------->  Remeshing at start of run step %-i. Remeshing iteration #%-i (#Ele=%i,#Nodes=%i) \n ',CtrlVar.CurrentRunStepNumber,JJ,MUAold.Nele,MUAold.Nnodes);
        end
        
        %  Determine new desired element sizes and identify elements for refinement
        %  or coarsening.
        [UserVar,RunInfo,Fold,xNod,yNod,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened]=...
            NewDesiredEleSizesAndElementsToRefineOrCoarsen2(UserVar,RunInfo,CtrlVar,MUAold,BCsOld,Fold,lold,Fold.GF,RuvOld,Lubvb);
        
        %  Remesh: either global-remeshing or local-mesh refinement.
        [UserVar,RunInfo,CtrlVar,MUAnew]=...
            Remeshing(UserVar,RunInfo,CtrlVar,MUAold,BCsOld,Fold,lold,Fold.GF,...
            xNod,yNod,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened);
        
        
        if MUAnew.Nele==0
            fprintf('No elements left in mesh! \n ')
            return
        end
        
        
        nNewElements=MUAnew.Nele-MUAold.Nele;
        nNewNodes=MUAnew.Nnodes-MUAold.Nnodes;
        
       
        [UserVar,RunInfo,Fnew,BCsNew,lnew]=MapFbetweenMeshes(UserVar,RunInfo,CtrlVar,MUAold,MUAnew,Fold,BCsOld,lold);
       
        
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
            hold on ;  [xGL,yGL]=PlotGroundingLines(CtrlVar,MUAold,Fold.GF,GLgeo,xGL,yGL,'r');
            title(sprintf('Adapt-mesh iteration #%i \n Before remeshing  \t #Ele=%-i, #Nodes=%-i, #nod=%-i',JJ,MUAold.Nele,MUAold.Nnodes,MUAold.nod))
            axis tight
            
            subplot(2,1,2)
            hold off
            xGL=[] ; yGL=[]; GLgeo=[];
            CtrlVar.PlotGLs=1;
            PlotMuaMesh(CtrlVar,MUAnew,[],CtrlVar.MeshColor);
            title(sprintf('After remeshing  \t #Ele=%-i, #Nodes=%-i, #nod=%-i',MUAnew.Nele,MUAnew.Nnodes,MUAnew.nod))
            hold on ;  [xGL,yGL]=PlotGroundingLines(CtrlVar,MUAnew,Fnew.GF,GLgeo,xGL,yGL,'r');
            axis tight
            
            
            
        end
        %%
        
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
                
    end
end




if ~isempty(CtrlVar.SaveAdaptMeshFileName)
    MUA=MUAnew;
    save(CtrlVar.SaveAdaptMeshFileName,'MUA') ;
    fprintf(CtrlVar.fidlog,'New mesh was saved in %s .\n',CtrlVar.SaveAdaptMeshFileName);
end

%%
%% map variables to new mesh

[UserVar,RunInfo,Fnew,BCsNew,lnew]=MapFbetweenMeshes(UserVar,RunInfo,CtrlVar,MUAold,MUAnew,Fold,BCsOld,lold);

%%

if CtrlVar.ManuallyDeactivateElements
    
    
    if CtrlVar.InfoLevelAdaptiveMeshing>=1
        fprintf('Manual deactivation of elements.\n'); 
    end
    
    
    xEle=Nodes2EleMean(MUAnew.connectivity,MUAnew.coordinates(:,1));
    yEle=Nodes2EleMean(MUAnew.connectivity,MUAnew.coordinates(:,2));
    ElementsToBeDeactivated=false(MUAnew.Nele,1);
    
    [UserVar,ElementsToBeDeactivated]=...
        DefineElementsToDeactivate(UserVar,RunInfo,CtrlVar,MUAnew,xEle,yEle,ElementsToBeDeactivated,Fnew.s,Fnew.b,Fnew.S,Fnew.B,Fnew.rho,Fnew.rhow,Fnew.ub,Fnew.vb,Fnew.ud,Fnew.vd,Fnew.GF);
    
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
    
   
    [UserVar,RunInfo,Fnew,BCsNew,lnew]=MapFbetweenMeshes(UserVar,RunInfo,CtrlVar,MUAold,MUAnew,Fold,BCsOld,lold);

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




%%
%save TestSave
%[UserVar,RunInfo,F,l,BCs,GF,dt]=uvh(UserVar,RunInfo,CtrlVar,MUAnew,Fnew,Fnew,lnew,lnew,BCsNew);
%%



end

