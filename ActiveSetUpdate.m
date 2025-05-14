



function [UserVar,RunInfo,BCs1,lambdahpos,isActiveSetModified,isActiveSetCyclical,Activated,DeActivated]=ActiveSetUpdate(UserVar,RunInfo,CtrlVar,MUA,F1,l1,BCs1,iActiveSetIteration,LastDeactivated,LastActivated)


narginchk(10,10)
nargoutchk(8,8)

RunInfo.Forward.ActiveSetConverged=1;
RunInfo.Forward.uvhIterationsTotal=0;
isActiveSetCyclical=NaN;

BCs1Input=BCs1;

if ~CtrlVar.ThicknessConstraints
    return
end


%  Thickness constraints used

%    NodesFixed: holds the nodal numbers nodes in the active set
%      ihactive: number nodes in active set
%
% -If min thickness is significantly greater than CtrlVar.ThickMin, eliminate all thickness constraints
% -If there are no thickness constraints, but some thicknesses (h1) are smaller than CtrlVar.ThickMin, introduce a new initial active set
%           and set velocity estimates to zero, and thickness to ThickMin at nodes in the active set
% Enter loop:
%      -Where h0 is at ThickMin set h1 to ThickMin and u1,v1 to u0, v0.  This modification is important for convergence
%      -Solve system by calling uvh2D
%      -If solution did not converge use an inner loop that tries to reset initial estimates for u1,v1 and h1.
%      -Based on sign of the Lagrange multipliers, take nodes out of active set
%      -Where h1<=ThickMin add those to active set, provided 1) they are not already in the active set and 2) they are to be taken out
%
%
% active set method is used to enforce thickness constraints


% new approach: just keep track of new hFixedNode
%
% have new fields in object BoundaryConditions
% hPosNode
% hPosValue



%%  Make sure that no thickness constraints have been added to nodes with user-defined boundary conditions

BCs1.hPosNode=setdiff(BCs1.hPosNode,BCs1.hFixedNode) ;

%%

if CtrlVar.ThicknessConstraintsInfoLevel>=1
    if numel(BCs1.hPosNode)>0 || CtrlVar.ThicknessConstraintsInfoLevel>=10
        fprintf(CtrlVar.fidlog,' Number of active thickness constraints is %-i \n',numel(BCs1.hPosNode));
    end
end

%%  Here we begin:  I have an active set already, now need to update it
% keep a copy of the old active set
LastActiveSet=BCs1.hPosNode;

%must now find the lambda values corresponding to nodes that where constrained to positive thickness

% hLambda=l1.h;
% lambdahpos=hLambda(numel(BCs1.hFixedNode)+numel(BCs1.hTiedNodeA)+1:end) ;%  I always put the hPos constraints at end of all other h constraints
lambdahpos=l1.h(numel(BCs1.hFixedNode)+numel(BCs1.hTiedNodeA)+1:end) ;%  I always put the hPos constraints at end of all other h constraints
% then this will work

%%  Mapping into 'physical' nodal basis if required
% remember lambdahpos=hLambda(numel(BCs1.hFixedNode)+numel(BCs1.hTiedNodeA)+1:end)
% actually most likely only need to do this if numel(hPosNode)>0
% If the L matrix was not in the FE basis, I simply calculate the reactions.
% The Reactions are calculated correctly irrespective of how the
%
if ~CtrlVar.LinFEbasis
    if numel(BCs1.hPosNode) >0

        Reactions=CalculateReactions(CtrlVar,MUA,BCs1,l1);
      
        ah=-Reactions.h./(F1.rho.*F1.dt) ;
        ah=ah(BCs1.hPosNode);
        lambdahpos=Reactions.h(BCs1.hPosNode);

    else
        lambdahpos=[];
        ah=[];
    end
end
%%


if numel(lambdahpos) ~= numel(BCs1.hPosNode)
    save TestSave ; error(' # of elements in lambdahpos must equal # of elements in BCs1.hPosNode')
end


%% Print information the solution for the Lagrange multipliers

if CtrlVar.ThicknessConstraintsInfoLevel>=10
    [~,I]=sort(lambdahpos);  % print out fixed nodes in the order of increasing lambda values
    fprintf(CtrlVar.fidlog,' Pos. thick. contraints: ')   ;
    fprintf(CtrlVar.fidlog,' \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \n \t \t \t \t \t \t',BCs1.hPosNode(I));
    fprintf(CtrlVar.fidlog,'\n   Lagrange multipliers: ') ;
    fprintf(CtrlVar.fidlog,' \t %+9.0g \t %+9.0g \t %+9.0g \t %+9.0g \t %+9.0g \t %+9.0g \t %+9.0g \t %+9.0g \t %+9.0f \t %+9.0f \n \t \t \t \t \t \t',lambdahpos(I)) ;


    fprintf(CtrlVar.fidlog,'\n            mass flux: ') ;

    fprintf(CtrlVar.fidlog,' \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \n \t \t \t \t \t \t',ah(I)) ;


    fprintf(CtrlVar.fidlog,'\n');
end



%% Updating the active set
% Now I've solved for h1 and if needed a new active set must be defined
%
% The new active set contains all nodes where h1 less than hmin that were not in the previous active set
% Those of the nodes in the previous active set with positive slack values
% Nodes in the previous set with negative slack values must be taken out of the set
% if the active-set method is selected, update active set
% The active set is created/modified and the problem solved again if the active set has changed


% Do I need to deactivate some thickness constraints?
% if any of the lambdahpos are positive, then these constraints must be deactivated

if numel(BCs1.hPosNode)>0   % are there any min thickness constraints? If so see if some should be deactivated

    % I divide here with rho and dt for this to have the same units as the mass balance
    % (distance/time)
    % 
    % isNegavtiveMassFluxSmall=lambdahpos>CtrlVar.ThicknessConstraintsLambdaPosThreshold./(F1.rho(BCs1.hPosNode).*F1.dt);  % if any of the Lagrange multipliers `lambdahpos' are positive, then these should be deactivated

    % Clearly only inactivate nodes if the mass flux needed to keep them active (ah) is negative.
    % But to also consider only inactivate if the negative flux is 
    % 
    %   ah < 0.01 hMin /dt 
    %
    % that is, if one were to stop subtracting this mass balance, then the thickness would increase to 0.01 above the minimum
    % thickness value over a time interval corresponding to one time unit.
    
    isNegavtiveMassFluxSmall=ah < -0.01*CtrlVar.ThickMin/F1.dt ; 

    NewInActiveConstraints=find(isNegavtiveMassFluxSmall);
    iNewInActiveConstraints=numel(NewInActiveConstraints);
    if iNewInActiveConstraints>0   % have any become inactive?
        BCs1.hPosNode(isNegavtiveMassFluxSmall)=[]; % Here I deactivate 

    end

else
    NewInActiveConstraints=[];
    iNewInActiveConstraints=numel(NewInActiveConstraints);
end

NodesDeactivated=LastActiveSet(NewInActiveConstraints);

% Do I need to activate some new thickness constraints?
%I=h1<=CtrlVar.ThickMin; % if thickness is less than ThickMin then further new thickness constraints must be introduced
I=F1.h<CtrlVar.ThickMin; % if thickness is less than ThickMin then further new thickness constraints must be introduced

NewActive=find(I);
NewActive=setdiff(NewActive,BCs1.hFixedNode) ;  % do not add active thickness constraints for nodes that already are included in the user-defined thickness boundary conditions,
                                                % even if this means that thicknesses at those nodes are less then MinThick
NewActive=setdiff(NewActive,BCs1.hPosNode);     % exclude those already in the active set
NewActive=setdiff(NewActive,NodesDeactivated);  % do not include those nodes at min thick that I now must release

if isfield(CtrlVar,"ActiveSet")
    if CtrlVar.ActiveSet.ExcludeNodesOfBoundaryElements
        BoundaryElementNodes=unique(MUA.connectivity([MUA.Boundary.Elements{:}]',:));
        NewActive=setdiff(NewActive,BoundaryElementNodes) ;
    end
end


iNewActiveConstraints=numel(NewActive);

if iNewActiveConstraints> CtrlVar.MaxNumberOfNewlyIntroducedActiveThicknessConstraints
    if CtrlVar.ThicknessConstraintsInfoLevel>=1
        fprintf(CtrlVar.fidlog,' Number of new active-set thickness constraints %-i larger then max number or newly added constraints %-i \n ',...
            iNewActiveConstraints,CtrlVar.MaxNumberOfNewlyIntroducedActiveThicknessConstraints);
        fprintf(CtrlVar.fidlog,' Only the smallest %-i thickness values are constrained \n',CtrlVar.MaxNumberOfNewlyIntroducedActiveThicknessConstraints);
    end
    [~,II]=sort(F1.h);
    NewActive=II(1:CtrlVar.MaxNumberOfNewlyIntroducedActiveThicknessConstraints);
    iNewActiveConstraints=CtrlVar.MaxNumberOfNewlyIntroducedActiveThicknessConstraints;

end

if iNewActiveConstraints>0
    if CtrlVar.ThicknessConstraintsInfoLevel>=1
        fprintf(CtrlVar.fidlog,' %i new active constraints \n',iNewActiveConstraints);
    end

    BCs1.hPosNode=[BCs1.hPosNode;NewActive] ; 
end


%% Consider adding LSF F1.LSFMask.NodesOut to thickness constraints

if CtrlVar.LevelSetMethod && CtrlVar.LevelSetMethodThicknessConstraints

    if isempty(F1.LSFMask)  %
        F1.LSFMask=CalcMeshMask(CtrlVar,MUA,F1.LSF,0);
    end

    LSFhPosNode=find(F1.LSFMask.NodesOut);
    LSFhAdditionalPosNodes=setdiff(LSFhPosNode,BCs1.hPosNode) ;

    if CtrlVar.ThicknessConstraintsInfoLevel>=1
        fprintf(' %i level-set thickness constraints \n',numel(LSFhAdditionalPosNodes))
    end

    BCs1.hPosNode=union(BCs1.hPosNode,LSFhPosNode);
    


end

if isfield(CtrlVar,"ActiveSet")
    if CtrlVar.ActiveSet.ExcludeNodesOfBoundaryElements
        BoundaryElementNodes=unique(MUA.connectivity([MUA.Boundary.Elements{:}]',:));
        BCs1.hPosNode=setdiff(BCs1.hPosNode,BoundaryElementNodes) ;
    end
end




% LastActiveSet is simply a copy of hPosNode at the beginning of the call
DeActivated=setdiff(LastActiveSet,BCs1.hPosNode)   ; % nodes in last active set that are no longer in the new one
Activated=setdiff(BCs1.hPosNode,LastActiveSet)  ; % nodes in new active set that were not in the previous one

nDeactivated=numel(DeActivated);
nActivated=numel(Activated);

% No further changes made to active set, EXEPT if the active set has become cyclical (see below) in which case the
% nodes labeled for deactivation are not deactivated

BCs1.hPosNodeDeActivated=DeActivated;  % I'm not really using this at the moment, but in the future it might be best to use this to determine if set has become cyclical 
BCs1.hPosNodeActivated=Activated;

%% print information on new active set
if CtrlVar.ThicknessConstraintsInfoLevel>=1
    if nDeactivated > 0 || nActivated > 0
        fprintf(CtrlVar.fidlog,'\n  Updating pos. thickness constraints: deactivated: %-i,  activated: %-i, total number of thickness constrains: %-i \n',...
            nDeactivated,nActivated,numel(BCs1.hPosNode));
        fprintf(CtrlVar.fidlog,'  Nodes inactivated: ')   ;
        fprintf(CtrlVar.fidlog,' \t %7i \t %7i \t %7i \t %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \n \t \t \t \t \t',DeActivated);

        fprintf(CtrlVar.fidlog,'\n    Nodes activated: ')   ;
        fprintf(CtrlVar.fidlog,' \t %7i \t %7i \t %7i \t %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \n \t \t \t \t \t',Activated);
        fprintf(CtrlVar.fidlog,'\n ')   ;
    else
        fprintf(CtrlVar.fidlog,'No pos.-thickness constraints activated or deactivated. \n')   ;
    end

end


%% check if set has become cyclical
%LastDeactivated=Deactivated; LastActivated=Activated;



isActiveSetCyclical=false ;



ActivatedAndPreviouslyDeactivatedDifference=setxor(Activated,LastDeactivated); % if empty then the sets of activated and previously de-activated nodes are identical
DeactivatedAndPreviouslyActivatedDifference=setxor(DeActivated,LastActivated); % if empty then the sets of de-activated and previously activated nodes are identical

if CtrlVar.ThicknessConstraintsInfoLevel>=1
    if ~isempty(LastDeactivated)
        if isempty(ActivatedAndPreviouslyDeactivatedDifference)
            fprintf(' Active-set: The set of nodes being activated is identical to the one previously deactivated. \n')
        end
    end

    if ~isempty(LastActivated)
        if isempty(DeactivatedAndPreviouslyActivatedDifference)
            fprintf(' Active-set: The set of nodes being deactivated is identical to the one previously activated. \n')

        end
    end
end

if ~(isempty(Activated) && isempty(DeActivated))  
    if isempty(ActivatedAndPreviouslyDeactivatedDifference)  && isempty(DeactivatedAndPreviouslyActivatedDifference)
        isActiveSetCyclical=true ;
        fprintf(' Active-set is cyclical. \n')
    end
end

% I now have a dilemma, since the set has become cyclical it is
% clear that if I deactivate any new nodes the thickness at
% those nodes will become too small in the next active-set
% interaction. A solution is simply not to deactivate and to add
% the cyclically deactivated nodes to the active set.

if isActiveSetCyclical
    BCs1.hPosNode=[BCs1.hPosNode;DeActivated] ;
end

%% Set the hPosValues, only need to do this once at the end
BCs1.hPosValue=BCs1.hPosNode*0+CtrlVar.ThickMin;
%%

if nDeactivated> 0 || nActivated>0
    if nDeactivated<CtrlVar.MinNumberOfNewlyIntroducedActiveThicknessConstraints && nActivated<CtrlVar.MinNumberOfNewlyIntroducedActiveThicknessConstraints
        fprintf("ActiveSetInitialisation: Not introducing any new thickness constraints as:\n")
        fprintf("\t #Deactivated=%i and #activated=%i nodes, both less than CtrlVar.MinNumberOfNewlyIntroducedActiveThicknessConstraints=%i. \n",nDeactivated,nActivated,CtrlVar.MinNumberOfNewlyIntroducedActiveThicknessConstraints)
        BCs1=BCs1Input;
        Activated=[];
        DeActivated=[];

    end
end




ChangeInActiveSet=setxor(BCs1.hPosNode,LastActiveSet) ;
nChangeInActiveSet=numel(ChangeInActiveSet);
%if nChangeInActiveSet == 0  ||  isActiveSetCyclical  % 1 Sept 2024: Decided that if active set changed, then
%isActiveSetModified=true, even if the active set has become cyclical. 
if nChangeInActiveSet == 0  
    isActiveSetModified=false;
    fprintf("ActiveSetUpdate: Active set not modified.\n")
else
    isActiveSetModified=true;
    fprintf("ActiveSetUpdate: Active set modified.\n")
end


% setxor(BCs1.hPosNode,BCs1Input.hPosNode)


RunInfo.Forward.uvhActiveSetIterations(CtrlVar.CurrentRunStepNumber)=iActiveSetIteration-1 ;
RunInfo.Forward.uvhActiveSetCyclical(CtrlVar.CurrentRunStepNumber)=isActiveSetCyclical;
RunInfo.Forward.uvhActiveSetConstraints(CtrlVar.CurrentRunStepNumber)=numel(BCs1.hPosNode);



end

