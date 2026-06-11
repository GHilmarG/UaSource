
function [UserVar,RunInfo,MUAnew]=MeshAdvanceRetreat(UserVar,RunInfo,CtrlVar,MUAold,BCsOld,Fold,lold,GFold,Ruv,Lubvb)


persistent MUA_Background

%%  FE mesh advance or retreat
if CtrlVar.InfoLevelAdaptiveMeshing>=1
    fprintf(CtrlVar.fidlog,' FE mesh advance or retreat, element activated/deactivated at time %-g \n ',CtrlVar.time);
end

if isempty(MUA_Background)
    
    if    exist(fullfile(cd,CtrlVar.FEmeshAdvanceRetreatBackgroundMeshFileName), 'file')  == 2 ...
            || exist(fullfile(cd,[CtrlVar.FEmeshAdvanceRetreatBackgroundMeshFileName,'.mat']), 'file')  == 2
        
        
        ListOfVariables=who('-file',CtrlVar.FEmeshAdvanceRetreatBackgroundMeshFileName);
        
        if CtrlVar.InfoLevelAdaptiveMeshing>=1
            fprintf('Reading ''background'' MUA from: %s ',CtrlVar.FEmeshAdvanceRetreatBackgroundMeshFileName)
        end
        
        if ismember('MUA_Background',ListOfVariables)
            
            load(CtrlVar.FEmeshAdvanceRetreatBackgroundMeshFileName,'MUA_Background')
            
        elseif ismember('MUA',ListOfVariables)
            
            
            VAR='MUA';
            MUA_Background=load(CtrlVar.FEmeshAdvanceRetreatBackgroundMeshFileName,VAR);
            MUA_Background=MUA_Background.(VAR);
            
        else
            
            fprintf(' Neither MUA or MUA_Background found in %s.! \n',CtrlVar.FEmeshAdvanceRetreatBackgroundMeshFileName)
            error('AdaptMesh:MUA_BackgroundNotFound','Where is MUA_Background?')
            
        end
        
        MUA_Background=UpdateMUA(CtrlVar,MUA_Background);
        if CtrlVar.InfoLevelAdaptiveMeshing>=1
            fprintf('done \n ')
        end
        
    else
        fprintf('File with background meshfile %s could not be found \n',CtrlVar.FEmeshAdvanceRetreatBackgroundMeshFileName)
        fprintf('Note: When using the''FEmeshAdvanceRetreat'' option, a background mesh must be defined.\n')
        fprintf('The background mesh is usually same as the intial mesh used in the run, but must be defined by the user.\n')
        fprintf('Possible solution: \n')
        fprintf('Load a MUA structure previously created in an initial meshing step and save in a seperate file. \n')
        fprintf('For example: load NewMeshFile.mat ; MUA_Background=MUA ; save(''Backgroundmesh'',''MUA_Background'') \n')
        fprintf('And in Ua2D_InitialUserInput: CtrlVar.FEmeshAdvanceRetreatBackgroundMeshFileName=''Backgroundmesh'' ; \n')
        fprintf('                              CtrlVar.FEmeshAdvanceRetreat=1 ; \n')
        
        error('Ua:AdaptMesh:NoBackgroundMeshfile',' exiting ' )
    end
    
end


hBackground=MapNodalVariablesFromMesh1ToMesh2(CtrlVar,MUAold,MUA_Background.coordinates(:,1),MUA_Background.coordinates(:,2),CtrlVar.ThickMin,Fold.h);
[UserVar,iDeactivatedElements]=FindElementsToDeactivate(UserVar,CtrlVar,MUA_Background,hBackground);


fprintf(CtrlVar.fidlog,'%i elements of background mesh deactivated. ',numel(find(iDeactivatedElements)));


if CtrlVar.doplots && CtrlVar.doAdaptMeshPlots && CtrlVar.InfoLevelAdaptiveMeshing>=10
    
    figure(10000)
    subplot(2,1,1,'replace')
    hold off
    PlotNodalBasedQuantities(MUA_Background.connectivity,MUA_Background.coordinates,hBackground,CtrlVar);
    xlabel('x (km)') ; ylabel('y (km)') ;
    colorbar ; title(colorbar,'(m)') ; % caxis([0 100])
    hold on
    
    title(sprintf('hBackground at t=%#5.1f ',CtrlVar.time))
    %tt=daspect ; daspect([mean(tt(1)+tt(2)) mean(tt(1)+tt(2)) tt(3)]);
    axis equal tight
    
    ind=hBackground<=CtrlVar.ThickMinDeactivateElements;
    plot(MUA_Background.coordinates(ind,1)/CtrlVar.PlotXYscale,MUA_Background.coordinates(ind,2)/CtrlVar.PlotXYscale,'.r')
    hold off
    
    
    subplot(2,1,2,'replace')
    hold off
    PlotFEmesh(MUA_Background.coordinates,MUA_Background.connectivity,CtrlVar)
    hold on
    ind=hBackground<=CtrlVar.ThickMinDeactivateElements;
    plot(MUA_Background.coordinates(ind,1)/CtrlVar.PlotXYscale,MUA_Background.coordinates(ind,2)/CtrlVar.PlotXYscale,'.r')
    
    xEleBackground=Nodes2EleMean(MUA_Background.connectivity,MUA_Background.coordinates(:,1));
    yEleBackground=Nodes2EleMean(MUA_Background.connectivity,MUA_Background.coordinates(:,2));
    
    plot(xEleBackground(iDeactivatedElements)/CtrlVar.PlotXYscale,yEleBackground(iDeactivatedElements)/CtrlVar.PlotXYscale,'og')
    title(sprintf('time=%-g \t #Ele=%-i, #Nodes=%-i, #nod=%-i, \n # EleBG=%-i, #NodesBG=%-i, #nodBG=%-i \n Green=Elements to deactivate',...
        CtrlVar.time,MUAold.Nele,MUAold.Nnodes,MUAold.nod,MUA_Background.Nele,MUA_Background.nod,MUA_Background.Nnodes))
    axis equal tight
end

% eliminate elements
if numel(find(iDeactivatedElements))>0
    %[coordinates,connectivity]=DeactivateElements(CtrlVar,iDeactivatedElements,MUA_Background.coordinates,MUA_Background.connectivity);
    MUAnew=DeactivateMUAelements(CtrlVar,MUA_Background,iDeactivatedElements);
    fprintf(CtrlVar.fidlog,'Change in the number of elements : %+i \n ',MUAnew.Nele-MUAold.Nele);
else
    MUAnew=MUAold;
end



end



