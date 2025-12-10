



function [UserVar,RunInfo,F1,l1,BCs1,isActiveSetModified,Activated,DeActivated]=ActiveSetInitialisation(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1)


narginchk(8,8)

 
BCs1Input=BCs1 ; 



%%  Make sure that no thickness constraints have been added to nodes with user-defined boundary conditions
% or because the user-defined boundary conditions have changed the last active set was updated/created.

BCs1.hPosNode=setdiff(BCs1.hPosNode,BCs1.hFixedNode) ;


%%

LastActiveSet=BCs1.hPosNode;

if CtrlVar.ThicknessConstraintsInfoLevel>=10
    fprintf(CtrlVar.fidlog,'  Enforcing min thickness of %-g using active-set method \n',CtrlVar.ThickMin);
end

%% Special case:  Check if all thicknesses are positive, in which case all thickness constraints should be deactivated
if min(F1.h) > 1.1*CtrlVar.ThickMin
    if CtrlVar.ThicknessConstraintsInfoLevel>=10
        fprintf(CtrlVar.fidlog,' Eliminating any possible previous thickness constraints as min(h1)=%-g>1.1*CtrlVar.ThickMin=%-g \n',min(F0.h),CtrlVar.ThickMin);
    end
    BCs1.hPosNode=[] ; 
end

%% Special case:  Check if there are no previous thickness constraints, in which case new should be introduced based on ice thickness
% possibly all thickness constraints were eliminated in AdapMesh when deactivating/activating elements
% so I check if there are no thickness constrains but h1 is at the min thick.
% if so then I introduce an initial active set based on h1
%!if isempty(Lhpos)

if isempty(BCs1.hPosNode)
    Active=find(F1.h<=CtrlVar.ThickMin);     % Although I might only want to add nodes to the active set if thickness is somewhat less than MinThick, I might here
                                             % have a situation where this is a restart run where the active set has been set to empty, or a run after re-meshing, and I
                                             % want the initial active sets to be similar to previous active set. 
                                             %
    Active=setdiff(Active,BCs1.hFixedNode) ; % do not add active thickness constraints for nodes that already are included in the user-defined thickness boundary conditions, 
                                             % even if this means that thicknesses at those nodes are less then MinThick
    BCs1.hPosNode=Active ;
    
    if numel(BCs1.hPosNode)>0
        if CtrlVar.ThicknessConstraintsInfoLevel>=1
            fprintf(CtrlVar.fidlog,' Introducing %-i new initial active constraints based on h1  \n', numel(BCs1.hPosNode));
        end
    end
end


% II= (F0.h<=CtrlVar.ThickMin) | (F1.h<=CtrlVar.ThickMin) ;
% F1.h(II)=CtrlVar.ThickMin;  F1.ub(II)=F0.ub(II) ; F1.vb(II)=F0.vb(II) ;  % modify initial guess for h1, POSSIBLY important for convergence
                                                                           % this is now done as ahead of uvh solve as iterate is made feasible

% However, I concluded (17 Dec, 2018) that it was better not to do this, as this can significantly increase the number of NR iterations.
% Better to do this only if the iteration does not converge.
% And then again on 17 Jan, 2019, it was found that it's better to keep this, as not doing so was found to have adverse effects on
% NR convergence rate. The reason for this is a bit unclear. This happened in after remeshing step and that may play a role. Anyhow,
% decided to revert back to previous tried-and-tested approach.




%% Consider adding LSF F1.LSFMask.NodesOut to thickness constraints

if CtrlVar.LevelSetMethod && CtrlVar.LevelSetMethodThicknessConstraints

    if isempty(F1.LSFMask)  % If I have already solved the LSF equation, this will not be empty and does not need to be recalculated (ToDo)
        F1.LSFMask=CalcMeshMask(CtrlVar,MUA,F1.LSF,0);
    end

    LSFhPosNode=find(F1.LSFMask.NodesOut);
    %LSFhAdditionalPosNodes=setdiff(BCs1.hPosNode,LSFhPosNode) ; 
    LSFhAdditionalPosNodes=setdiff(LSFhPosNode,BCs1.hPosNode) ; 

   if CtrlVar.ThicknessConstraintsInfoLevel>=1
        fprintf(' %i new level set active thickness constraints introduced \n',numel(LSFhAdditionalPosNodes))
    end

    BCs1.hPosNode=union(BCs1.hPosNode,LSFhPosNode); 
    
    
end


% I think the way LastActiveSet is initialized, both Deactivated  will always be empty by construction 
DeActivated=setdiff(LastActiveSet,BCs1.hPosNode)   ; % nodes in last active set that are no longer in the new one
Activated=setdiff(BCs1.hPosNode,LastActiveSet)  ; % nodes in new active set that were not in the previous one

nDeactivated=numel(DeActivated);
nActivated=numel(Activated);

%%


if nDeactivated> 0 || nActivated>0
    if nDeactivated<CtrlVar.MinNumberOfNewlyIntroducedActiveThicknessConstraints && nActivated<CtrlVar.MinNumberOfNewlyIntroducedActiveThicknessConstraints
        fprintf("ActiveSetInitialisation: Not introducing any new thickness constraints as:\n")
        fprintf("\t #Deactivated=%i and #activated=%i nodes, both less than CtrlVar.MinNumberOfNewlyIntroducedActiveThicknessConstraints=%i. \n",nDeactivated,nActivated,CtrlVar.MinNumberOfNewlyIntroducedActiveThicknessConstraints)

        BCs1=BCs1Input;
   
        Activated=[];
        DeActivated=[];
        nDeactivated=0;
        nActivated=0;
    end
end

%%

if nDeactivated==0   && nActivated==0
    isActiveSetModified=false;
    fprintf("ActiveSetInitialisation: Active set not modified.\n")
else
    isActiveSetModified=true;
    fprintf("ActiveSetInitialisation: Active set modified.\n")
end


%%
BCs1.hPosNodeDeActivated=DeActivated;  % I'm not really using this at the moment, but in the future it might be best to use this to determine if set has become cyclical 
BCs1.hPosNodeActivated=Activated;

%% Set the hPosValues, only need to do this once at the end
BCs1.hPosValue=BCs1.hPosNode*0+CtrlVar.ThickMin;
%%




end