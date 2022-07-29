



function [UserVar,RunInfo,BCs1,lambdahpos,isActiveSetModified,isActiveSetCyclical,Activated,Released]=ActiveSetUpdate(UserVar,RunInfo,CtrlVar,MUA,F1,l1,BCs1,iActiveSetIteration,LastReleased,LastActivated)


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

%    NodesFixed: holdes the nodal numbers nodes in the active set
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



%%




if CtrlVar.ThicknessConstraintsInfoLevel>=1
    if numel(BCs1.hPosNode)>0 || CtrlVar.ThicknessConstraintsInfoLevel>=10
        fprintf(CtrlVar.fidlog,' Number of active thickness constraints is %-i \n',numel(BCs1.hPosNode));
    end
end

%%  Here we begin:  I have an active set already, now need to update it
% keep a copy of the old active set
LastActiveSet=BCs1.hPosNode;

%must now find the lambda values corresponding to nodes that where constrainted to pos thickness

% hLambda=l1.h;
% lambdahpos=hLambda(numel(BCs1.hFixedNode)+numel(BCs1.hTiedNodeA)+1:end) ;%  I always put the hPos constraints at end of all other h constraints
lambdahpos=l1.h(numel(BCs1.hFixedNode)+numel(BCs1.hTiedNodeA)+1:end) ;%  I always put the hPos constraints at end of all other h constraints
% then this will work

%%  Mapping into 'physical' nodal basis if required
% remember lambdahpos=hLambda(numel(BCs1.hFixedNode)+numel(BCs1.hTiedNodeA)+1:end)
% actually most likely only need to do this if numel(hPosNode)>0
% If the L matrix was not in the FE basis, I simply calculate the reactions.
% The Reactions are calculated correctly irrespectivly of how the
%
if ~CtrlVar.LinFEbasis
    if numel(BCs1.hPosNode) >0

        Reactions=CalculateReactions(CtrlVar,MUA,BCs1,l1);
        lambdahpos=Reactions.h(BCs1.hPosNode);

    else
        lambdahpos=[];
    end
end
%%


if numel(lambdahpos) ~= numel(BCs1.hPosNode)
    save TestSave ; error(' # of elements in lambdahpos must equal # of elements in BCs1.hPosNode')
end


%% Print information the solution for the Lagrange multipliers

if CtrlVar.ThicknessConstraintsInfoLevel>=10
    [~,I]=sort(lambdahpos);  % print out fixed nodes in the order of increasing lambda values
    fprintf(CtrlVar.fidlog,'            Nodes fixed: ')   ;
    fprintf(CtrlVar.fidlog,' \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \n \t \t \t \t \t \t',BCs1.hPosNode(I));
    fprintf(CtrlVar.fidlog,'\n   Lagrange multipliers: ') ;
    fprintf(CtrlVar.fidlog,' \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \n \t \t \t \t \t \t',lambdahpos(I)) ;
    fprintf(CtrlVar.fidlog,'\n');
end



%% Updating the active set
% Now I've solved for h1 and if needed a new active set must be defined
%
% The new active set contains all nodes where h1 less than hmin that were not in the previous active set
% Those of the nodes in the previous active set with positve slack values
% Nodes in the previous set with negative slack values must be taken out of the set
% if the active-set method is selected, update active set
% The active set is created/modified and the problem solved again if the active set has changed


% Do I need to in-activate some thickness constraints?
% if any of the lambdahpos are positive, then these constraints must be in-activated

if numel(BCs1.hPosNode)>0   % are there any thickness constraints? If so see if some should be in-activated

    % I divide here with rho for this to have the same units as the mass balance
    % (distance/time)
    I=lambdahpos>CtrlVar.ThicknessConstraintsLambdaPosThreshold./F1.rho(BCs1.hPosNode);  % if any of the Lagrange multipliers `lambdahpos' are positive, then these should be in-activated
    NewInActiveConstraints=find(I);
    iNewInActiveConstraints=numel(NewInActiveConstraints);
    if iNewInActiveConstraints>0   % have any become inactive?
        BCs1.hPosNode(I)=[]; BCs1.hPosValue(I)=[];
     
    end

else
    NewInActiveConstraints=[];
    iNewInActiveConstraints=numel(NewInActiveConstraints);
end

NodesReleased=LastActiveSet(NewInActiveConstraints);

% Do I need to activate some new thickness constraints?
%I=h1<=CtrlVar.ThickMin; % if thickness is less than ThickMin then further new thickness constraints must be introduced
I=F1.h<=(CtrlVar.ThickMin-100*eps); % if thickness is less than ThickMin then further new thickness constraints must be introduced

NodesWithTooSmallThick=find(I);
NewActive=setdiff(NodesWithTooSmallThick,BCs1.hPosNode);  % exclude those already in the active set
NewActive=setdiff(NewActive,NodesReleased);  % do not include those nodes at min thick that I now must release


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
   
    BCs1.hPosNode=[BCs1.hPosNode;NewActive] ; BCs1.hPosValue=BCs1.hPosNode*0+CtrlVar.ThickMin;
end


%% Consider adding LSF F1.LSFMask.NodesOut to thickness constraints

if CtrlVar.LevelSetMethod && CtrlVar.LevelSetMethodThicknessConstraints

    if isempty(F1.LSFMask)  % 
        F1.LSFMask=CalcMeshMask(CtrlVar,MUA,F1.LSF,0);
    end

    LSFhPosNode=find(F1.LSFMask.NodesOut);
    LSFhAdditionalPosNodes=setdiff(LSFhPosNode,BCs1.hPosNode) ; 

   if CtrlVar.ThicknessConstraintsInfoLevel>=1
        fprintf(' %i additional LSF active constraints \n',numel(LSFhAdditionalPosNodes))
    end

    BCs1.hPosNode=union(BCs1.hPosNode,LSFhPosNode); 
    BCs1.hPosValue=BCs1.hPosNode*0+CtrlVar.ThickMin;


end

Released=setdiff(LastActiveSet,BCs1.hPosNode)   ; % nodes in last active set that are no longer in the new one
Activated=setdiff(BCs1.hPosNode,LastActiveSet)  ; % nodes in new active set that were not in the previous one

nReleased=numel(Released);
nActivated=numel(Activated);

% modify initial guess for h1, important for convergence
%h1(NewActive)=ThickMin;


% %% print information on new active set
% if CtrlVar.ThicknessConstraintsInfoLevel>=1
%     if iNewInActiveConstraints> 0 || iNewActiveConstraints> 0
%         fprintf(CtrlVar.fidlog,' Updating pos. thickness constraints: in-activated: %-i,  activated: %-i, total number of thickness constrains: %-i \n',...
%             iNewInActiveConstraints,iNewActiveConstraints,numel(BCs1.hPosNode));
%         fprintf(CtrlVar.fidlog,'  Nodes inactivated: ')   ;
%         fprintf(CtrlVar.fidlog,' \t %7i \t %7i \t %7i \t %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \n \t \t \t \t \t',NodesReleased);
% 
%         fprintf(CtrlVar.fidlog,'\n    Nodes activated: ')   ;
%         fprintf(CtrlVar.fidlog,' \t %7i \t %7i \t %7i \t %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \n \t \t \t \t \t',NewActive);
%         fprintf(CtrlVar.fidlog,'\n ')   ;
%     else
%         fprintf(CtrlVar.fidlog,'No pos.-thickness constraints activated or deactivated. \n')   ;
%     end
% 
% end

%% print information on new active set
if CtrlVar.ThicknessConstraintsInfoLevel>=1
    if nReleased > 0 || nActivated > 0
        fprintf(CtrlVar.fidlog,'\n  Updating pos. thickness constraints: in-activated: %-i,  activated: %-i, total number of thickness constrains: %-i \n',...
            nReleased,nActivated,numel(BCs1.hPosNode));
        fprintf(CtrlVar.fidlog,'  Nodes inactivated: ')   ;
        fprintf(CtrlVar.fidlog,' \t %7i \t %7i \t %7i \t %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \n \t \t \t \t \t',Released);

        fprintf(CtrlVar.fidlog,'\n    Nodes activated: ')   ;
        fprintf(CtrlVar.fidlog,' \t %7i \t %7i \t %7i \t %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \n \t \t \t \t \t',Activated);
        fprintf(CtrlVar.fidlog,'\n ')   ;
    else
        fprintf(CtrlVar.fidlog,'No pos.-thickness constraints activated or deactivated. \n')   ;
    end

end


%% check if set has become cyclical
%LastReleased=Released; LastActivated=Activated;



isActiveSetCyclical=false ;



ActivatedAndPreviouslyReleasedDifference=setxor(Activated,LastReleased); % if empty then the sets of activated and previously de-activated nodes are identical
ReleasedAndPreviouslyActivatedDifference=setxor(Released,LastActivated); % if empty then the sets of de-activated and previously activated nodes are identical

if CtrlVar.ThicknessConstraintsInfoLevel>=1
    if ~isempty(LastReleased)
        if isempty(ActivatedAndPreviouslyReleasedDifference)
            fprintf(' Active-set: The set of nodes being activated is identical to the one previously in-activated. \n')
        end
    end

    if ~isempty(LastActivated)
        if isempty(ReleasedAndPreviouslyActivatedDifference)
            fprintf(' Active-set: The set of nodes being in-activated is identical to the one previously activated. \n')

        end
    end
end

if isempty(ActivatedAndPreviouslyReleasedDifference)  && isempty(ReleasedAndPreviouslyActivatedDifference)
    isActiveSetCyclical=true ;
    fprintf(' Active-set is cyclical. \n')
end

% I now have a dilemma, since the set has become cyclical it is
% clear that if I deactivate any new nodes the thickness at
% those nodes will become too small in the next active-set
% interation. A solution is simply not to deactivate and to add
% the cyclically deactivated nodes to the active set.

if isActiveSetCyclical
    BCs1.hPosNode=[BCs1.hPosNode;Released] ; 
    BCs1.hPosValue=BCs1.hPosNode*0+CtrlVar.ThickMin;
end




if nReleased> 0 || nActivated>0
    if nReleased<CtrlVar.MinNumberOfNewlyIntroducedActiveThicknessConstraints && nActivated<CtrlVar.MinNumberOfNewlyIntroducedActiveThicknessConstraints
        fprintf("ActiveSetInitialisation: Not introducing any new thickness constraints as:\n")
        fprintf("\t #released=%i and #activated=%i nodes, both less than CtrlVar.MinNumberOfNewlyIntroducedActiveThicknessConstraints=%i. \n",nReleased,nActivated,CtrlVar.MinNumberOfNewlyIntroducedActiveThicknessConstraints)
        BCs1=BCs1Input;
        Activated=[];
        Released=[];

    end
end




ChangeInActiveSet=setxor(BCs1.hPosNode,LastActiveSet) ;
nChangeInActiveSet=numel(ChangeInActiveSet);
if nChangeInActiveSet == 0  ||  isActiveSetCyclical
    isActiveSetModified=false;
    fprintf("ActiveSetUpdate: Active set not modified.\n")
else
    isActiveSetModified=true;
    fprintf("ActiveSetUpdate: Active set modified.\n")
end


RunInfo.Forward.uvhActiveSetIterations(CtrlVar.CurrentRunStepNumber)=iActiveSetIteration-1 ;
RunInfo.Forward.uvhActiveSetCyclical(CtrlVar.CurrentRunStepNumber)=isActiveSetCyclical;
RunInfo.Forward.uvhActiveSetConstraints(CtrlVar.CurrentRunStepNumber)=numel(BCs1.hPosNode);

