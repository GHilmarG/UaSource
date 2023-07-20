



function [UserVar,RunInfo,F1,l1,BCs1,isActiveSetModified,Activated,Released]=ActiveSetInitialisation(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1)


 
 BCs1Input=BCs1 ; 

LastActiveSet=BCs1.hPosNode;

if CtrlVar.ThicknessConstraintsInfoLevel>=10
    fprintf(CtrlVar.fidlog,'  Enforcing min thickness of %-g using active-set method \n',CtrlVar.ThickMin);
end

%% Special case:  Check if all thicknesses are positive, in which case all thickness constraints should be deactiated
if min(F1.h) > 1.1*CtrlVar.ThickMin
    if CtrlVar.ThicknessConstraintsInfoLevel>=10
        fprintf(CtrlVar.fidlog,' Eliminating any possible previous thickness constraints as min(h1)=%-g>1.1*CtrlVar.ThickMin=%-g \n',min(F0.h),CtrlVar.ThickMin);
    end
    BCs1.hPosNode=[] ; BCs1.hPosValue=[];
end

%% Special case:  Check if there are no previous thickness constraints, in which case new should be introduced based on ice thickness
% possibly all thickness constraints were eliminated in AdapMesh when deactivating/activating elements
% so I check if there are no thickness constrains but h0 is at the min thick.
% if so then I introduce an initial active set based on h0
%!if isempty(Lhpos)

if isempty(BCs1.hPosNode)
    Active=find(F1.h<=CtrlVar.ThickMin);
    BCs1.hPosNode=Active ; BCs1.hPosValue=BCs1.hPosNode*0+CtrlVar.ThickMin;
    % F1.ub(BCs1.hPosNode)=0 ; F1.vb(BCs1.hPosNode)=0; 
    F1.h(BCs1.hPosNode)=CtrlVar.ThickMin; % 
    if numel(BCs1.hPosNode)>0
        if CtrlVar.ThicknessConstraintsInfoLevel>=1
            fprintf(CtrlVar.fidlog,' Introducing %-i new initial active constraints based on h0  \n', numel(BCs1.hPosNode));
        end
    end
end


II= (F0.h<=CtrlVar.ThickMin) | (F1.h<=CtrlVar.ThickMin) ;

F1.h(II)=CtrlVar.ThickMin;  F1.ub(II)=F0.ub(II) ; F1.vb(II)=F0.vb(II) ;  % modify initial guess for h1, POSSIBLY important for convergence

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
        fprintf(' %i new LSF active constraints \n',numel(LSFhAdditionalPosNodes))
    end

    BCs1.hPosNode=union(BCs1.hPosNode,LSFhPosNode); 
    BCs1.hPosValue=BCs1.hPosNode*0+CtrlVar.ThickMin;
    
end



Released=setdiff(LastActiveSet,BCs1.hPosNode)   ; % nodes in last active set that are no longer in the new one
Activated=setdiff(BCs1.hPosNode,LastActiveSet)  ; % nodes in new active set that were not in the previous one

nReleased=numel(Released);
nActivated=numel(Activated);

%%


if nReleased> 0 || nActivated>0
    if nReleased<CtrlVar.MinNumberOfNewlyIntroducedActiveThicknessConstraints && nActivated<CtrlVar.MinNumberOfNewlyIntroducedActiveThicknessConstraints
        fprintf("ActiveSetInitialisation: Not introducing any new thickness constraints as:\n")
        fprintf("\t #released=%i and #activated=%i nodes, both less than CtrlVar.MinNumberOfNewlyIntroducedActiveThicknessConstraints=%i. \n",nReleased,nActivated,CtrlVar.MinNumberOfNewlyIntroducedActiveThicknessConstraints)

        BCs1=BCs1Input;
   
        Activated=[];
        Released=[];
        nReleased=0;
        nActivated=0;
    end
end

%%

if nReleased==0   && nActivated==0
    isActiveSetModified=false;
    fprintf("ActiveSetInitialisation: Active set not modified.\n")
else
    isActiveSetModified=true;
    fprintf("ActiveSetInitialisation: Active set modified.\n")
end








end