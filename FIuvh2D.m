function [UserVar,RunInfo,F1,l1,BCs1,dt]=FIuvh2D(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1)

% save TestSave
% load TestSave ; CtrlVar.LinFEbasis=false ; [UserVar,RunInfo,F1,l1,BCs1,dt]=FIuvh2D(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1);


narginchk(9,9)
nargoutchk(6,6)

dt=CtrlVar.dt;

RunInfo.Forward.ActiveSetConverged=1;
RunInfo.Forward.IterationsTotal=0;

if ~CtrlVar.ThicknessConstraints
    
    
    [UserVar,RunInfo,F1,l1,BCs1]=uvh2D(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);
    
    CtrlVar.NumberOfActiveThicknessConstraints=0;
    
else
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
    
    
    ActiveSetReset=0; VariablesReset=0;  % keep this out here, only set once, whereas ReduceTimeStep is reset in each active-set iteration.
    iCounter=1;
    %hfixednodeposNew=[];
    if CtrlVar.ThicknessConstraintsInfoLevel>=10
        fprintf(CtrlVar.fidlog,'  Enforcing min thickness of %-g using active-set method \n',CtrlVar.ThickMin);
    end
    
    if min(F1.h) > 1.1*CtrlVar.ThickMin
        if CtrlVar.ThicknessConstraintsInfoLevel>=10
            fprintf(CtrlVar.fidlog,' Eliminating any possible previous thickness constraints as min(h1)=%-g>1.1*CtrlVar.ThickMin=%-g \n',min(F0.h),CtrlVar.ThickMin);
        end
        BCs1.hPosNode=[] ; BCs1.hPosValue=[];
    end
    
    % possibly all thickness constraints were eliminated in AdapMesh when deactivating/activating elements
    % so I check if there are no thickness constrains but h0 is at the min thick.
    % if so then I introduce an initial active set based on h0
    %!if isempty(Lhpos)
    if isempty(BCs1.hPosNode)
        Active=find(F0.h<=CtrlVar.ThickMin);
        BCs1.hPosNode=Active ; BCs1.hPosValue=BCs1.hPosNode*0+CtrlVar.ThickMin;
        F1.ub(BCs1.hPosNode)=0 ; F1.vb(BCs1.hPosNode)=0; F1.h(BCs1.hPosNode)=CtrlVar.ThickMin; % now set estimates of velocities at hPosNode to zero
        if numel(BCs1.hPosNode)>0
            if CtrlVar.ThicknessConstraintsInfoLevel>=1
                fprintf(CtrlVar.fidlog,' Introducing %-i new initial active constraints based on h0  \n', numel(BCs1.hPosNode));
            end
        end
    end
    
    isActiveSetModified=1;
    it=0;
    
    nlIt=zeros(CtrlVar.ThicknessConstraintsItMax+1,1)*2+NaN;  %
    Released=[] ; Activated=[];
    iCountCyclicActiveSet=0;  % counts the number the same set of nodes is activated and then deactivated
    
    while true   % active-set loop
        
        if ~isActiveSetModified
            fprintf(' Leaving active-set loop because active set did not change in last active-set iteration. \n')
            break
            
        end
        
        
        it=it+1;
        isActiveSetModified=0;
        
        fprintf(CtrlVar.fidlog,' ----------> Active-set Iteration #%-i <----------  \n',it);
        
        if CtrlVar.ThicknessConstraintsInfoLevel>=1
            if numel(BCs1.hPosNode)>0 || CtrlVar.ThicknessConstraintsInfoLevel>=10
                fprintf(CtrlVar.fidlog,' Number of active thickness constraints is %-i \n',numel(BCs1.hPosNode));
            end
        end
        
        II= (F0.h<=CtrlVar.ThickMin) | (F1.h<=CtrlVar.ThickMin) ;
        F1.h(II)=CtrlVar.ThickMin;  F1.ub(II)=F0.ub(II) ; F1.vb(II)=F0.vb(II) ;  % modify initial guess for h1, POSSIBLY important for convergence
        % However, I concluded (17 Dec, 2018) that it was better not to do this, as this can significantly inctrease the number of NR iterations.
        % Better to do this only if the iteration does not converge.
        % And then again on 17 Jan, 2019, it was found that it's better to keep this, as not doing so was found to have adverse effects on
        % NR convergence rate. The reason for this is a bit unclear. This happened in after remeshing step and that may play a role. Anyhow,
        % decided to revert back to previous tried-and-tested approach.
        
        F1.h(BCs1.hPosNode)=CtrlVar.ThickMin;
        
        
        uvhIt=0; ReduceTimeStep=0;% needed in inner loop
        
        %% Solve uv for the current active set
        while uvhIt<=2  % I have a possible inner iteration here because
            
            uvhIt=uvhIt+1;
            % if uvh2D does not converge fully, it is usually best to just update the
            % active set on the partially converged solution
            
            %if any(h0<CtrlVar.ThickMin) ;
            %    save TestSave ; warning('h0 negative on input to uvh2D. min(h0)=%g ') ;
            %    keyboard
            %end
            
            
            
            if CtrlVar.DebugMode
                warning('FIuvh2D:nonconvergent','Ahead of uvh2D call')
                filename='Dump_Before_FIuvh2D.mat';
                fprintf('Saving all data in %s \n',filename)
                save(filename)
            end
            
            
            [UserVar,RunInfo,F1,l1,BCs1]=uvh2D(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);
            nlIt(iCounter)=RunInfo.Forward.Iterations;  iCounter=iCounter+1;
            RunInfo.Forward.Iterations=mean(nlIt,'omitnan');
            
            % keep a copy of the old active set
            LastActiveSet=BCs1.hPosNode;
            
            %must now find the lambda values corresponding to nodes that where constrainted to pos thickness
            
            hLambda=l1.h;
            %lambdahpos=hLambda(numel(BCs1.hFixedNode)+1:end);  % if I always put the hPos constraints at end of all other h constraints
            lambdahpos=hLambda(numel(BCs1.hFixedNode)+numel(BCs1.hTiedNodeA)+1:end) ;%  I always put the hPos constraints at end of all other h constraints
            % then this will work
            
            %%  Mapping ino 'physical' nodal basis if required
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
            
            if RunInfo.Forward.Converged==1
                break
            end
            
            if CtrlVar.DebugMode
                warning('FIuvh2D:nonconvergent','uvh2D did not converge')
                filename='Dump_FIuvh2D.mat';
                fprintf('Saving all data in %s \n',filename)
                save(filename)
                error('asdf')
            end
            
            % If not converged try:
            % 1) Reset variables to values at the beginning of time ste
            % 2) If that does not work, reset active set
            % 3) If that does not work reset time step.
            
            if ~VariablesReset
                
                VariablesReset=1;
                if CtrlVar.ThicknessConstraintsInfoLevel>=1
                    fprintf(CtrlVar.fidlog,' uvh solution did not converge. Creating a new uvh starting point for the solver. \n');
                end
                
                %                 F1.ub=F0.ub ; F1.vb=F0.vb; F1.ud=F0.ud ; F1.vd=F0.vd ; F1.h=F0.h;
                %                 F1.h(BCs1.hPosNode)=CtrlVar.ThickMin;  % consider adding
                %                 l1.ubvb=l1.ubvb*0 ; l1.udvd=l1.udvd*0; l1.h=l1.h*0;
                
                dtOld=dt;
                dt=dt/2; CtrlVar.dt=dt;
                fprintf(CtrlVar.fidlog,' Warning : Reducing time step from %-g to %-g \n',dtOld,CtrlVar.dt);
                
                F0.h(F0.h<=CtrlVar.ThickMin)=CtrlVar.ThickMin ; F0.h(BCs1.hPosNode)=CtrlVar.ThickMin;
                [UserVar,RunInfo,F0,l0]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs1,F0,l0);
                %[UserVar,dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F0,BCs1);
                [UserVar,dhdt]=dhdtExplicitSUPG(UserVar,CtrlVar,MUA,F0,BCs1);
                F1=F0;
                F1.h=F0.h+dhdt.*CtrlVar.dt ;
                F1.h(F1.h<=CtrlVar.ThickMin)=CtrlVar.ThickMin ; F1.h(BCs1.hPosNode)=CtrlVar.ThickMin;
                [F1.b,F1.s,F1.h,F1.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F1.h,F1.S,F1.B,F1.rho,F1.rhow);
                l1.ubvb=l1.ubvb*0 ; l1.udvd=l1.udvd*0; l1.h=l1.h*0;
                [UserVar,RunInfo,F1,l1]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs1,F1,l1);
                
                
                
            elseif ~ReduceTimeStep
                
                ReduceTimeStep=1;
                dtOld=CtrlVar.dt;
                dt=dt/2; CtrlVar.dt=dt;
                
                fprintf(CtrlVar.fidlog,' Warning : Reducing time step from %-g to %-g \n',dtOld,CtrlVar.dt);
                fprintf(CtrlVar.fidlog,'Also resetting u1, v1, h1 to ub0, vb0 and h0, and setting estimates for Lagrange parameters to zero. \n');
                
                
                if CtrlVar.WriteRunInfoFile
                    fprintf(RunInfo.File.fid,' Warning : Reducing time step from %-g to %-g \n',dtOld,CtrlVar.dt);
                    fprintf(RunInfo.File.fid,'Also resetting u1, v1, h1 to ub0, vb0 and h0, and setting estimates for Lagrange parameters to zero. \n');
                    
                end
                
                F1.ub=F0.ub ; F1.vb=F0.vb ;  F1.ud=F0.ud ; F1.vd=F0.vd ; F1.h=F0.h;
                l1.ubvb=l1.ubvb*0 ; l1.udvd=l1.udvd*0; l1.h=l1.h*0;
                
            end
        end  % uvh loop
        
        %% Print information the solution for the Lagrange multipliers
        
        if CtrlVar.ThicknessConstraintsInfoLevel>=10
            [~,I]=sort(lambdahpos);  % print out fixed nodes in the order of increasing lambda values
            fprintf(CtrlVar.fidlog,'            Nodes fixed: ')   ;
            fprintf(CtrlVar.fidlog,' \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \n \t \t \t \t \t \t',BCs1.hPosNode(I));
            fprintf(CtrlVar.fidlog,'\n   Lagrange multipliers: ') ;
            fprintf(CtrlVar.fidlog,' \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \n \t \t \t \t \t \t',hLambda(I)) ;
            fprintf(CtrlVar.fidlog,'\n');
        end
        
        %% Check exit criteria that must be done after uv calculation but before new update of the active set
        % On exit I want the active-set to be the one used when calculating the uv solution.
        %
        % Either: 1) the active-set has not been changed in which case I can exist at the end of the loop or the beginning of the new one,
        %     or  2) the active-set has changed but I still want to exit the loop due to some other criteria being met.
        %
        % For case 2, I must exit the loop here because I've already have calculated the
        % velocities, but I have not updated the active set.
        
        if  it > CtrlVar.ThicknessConstraintsItMax
            fprintf(' Leaving active-set pos. thickness loop because number of active-set iteration (%i) greater than maximum allowed (CtrlVar.ThicknessConstraintsItMax=%i). \n ',it,CtrlVar.ThicknessConstraintsItMax)
            break
        end
        
        if iCountCyclicActiveSet>=CtrlVar.ThicknessConstraintsItMaxCycles
            
            fprintf(' Leaving active-set pos. thickness loop because it has become cyclical with number of cycles (%i) equal to maximum allowed (CtrlVar.ThicknessConstraintsItMaxCycles=%i). \n ',iCountCyclicActiveSet,CtrlVar.ThicknessConstraintsItMaxCycles)
            break
            
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
            
            % sometimes constraints are being activated and in-activated over and over again. A possible remedy is not to in-activate constraints
            % immediately and to introduce a 1% threshold value
            
            lambdahposThreshold=CtrlVar.ThicknessConstraintsLambdaPosThreshold;
            %
            %             if it>4 && CtrlVar.ThicknessConstraintsLambdaPosThreshold==0
            %
            %                 lambdahposThreshold=-mean(lambdahpos(lambdahpos<0))/100;
            %                 if isnan(lambdahposThreshold) ; lambdahposThreshold=0 ; end
            %                 if CtrlVar.ThicknessConstraintsInfoLevel>=1 ;
            %                     fprintf(CtrlVar.fidlog,' Introducing a min threshold of  %-g for  in-activating thickness constraints. \n',lambdahposThreshold);
            %                 end
            %             end
            
            
            %%
            I=lambdahpos>lambdahposThreshold  ;  % if any of the Lagrange multipliers `lambdahpos' are positive, then these should be in-activated
            NewInActiveConstraints=find(I);
            iNewInActiveConstraints=numel(NewInActiveConstraints);
            if iNewInActiveConstraints>0   % have any become inactive?
                BCs1.hPosNode(I)=[]; BCs1.hPosValue(I)=[];
                isActiveSetModified=1;
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
            isActiveSetModified=1;
            BCs1.hPosNode=[BCs1.hPosNode;NewActive] ; BCs1.hPosValue=BCs1.hPosNode*0+CtrlVar.ThickMin;
        end
        
        
        % modify initial guess for h1, important for convergence
        %h1(NewActive)=ThickMin;
        
        
        %% print information on new active set
        if CtrlVar.ThicknessConstraintsInfoLevel>=1
            if iNewInActiveConstraints> 0 || iNewActiveConstraints> 0
                fprintf(CtrlVar.fidlog,' Updating pos. thickness constraints: in-activated: %-i,  activated: %-i, total number of thickness constrains: %-i \n',...
                    iNewInActiveConstraints,iNewActiveConstraints,numel(BCs1.hPosNode));
                fprintf(CtrlVar.fidlog,'  Nodes inactivated: ')   ;
                fprintf(CtrlVar.fidlog,' \t %7i \t %7i \t %7i \t %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \n \t \t \t \t \t',NodesReleased);
                
                fprintf(CtrlVar.fidlog,'\n    Nodes activated: ')   ;
                fprintf(CtrlVar.fidlog,' \t %7i \t %7i \t %7i \t %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \n \t \t \t \t \t',NewActive);
                fprintf(CtrlVar.fidlog,'\n ')   ;
            else
                fprintf(CtrlVar.fidlog,'No pos.-thickness constraints activated or deactivated. \n')   ;
            end
            
            if isActiveSetModified==1 && it <= CtrlVar.ThicknessConstraintsItMax
                fprintf(' Active set modified.\n')
            end
        end
        
        %% check if set has become cyclical
        LastReleased=Released; LastActivated=Activated;
        Released=setdiff(LastActiveSet,BCs1.hPosNode)   ;% nodes in last active set that are no longer in the new one
        Activated=setdiff(BCs1.hPosNode,LastActiveSet)  ;% nodes in new active set that were not in the previous one
        
        isActiveSetCyclical=false ;
        
        if it>1
            
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
                iCountCyclicActiveSet=iCountCyclicActiveSet+1;
                isActiveSetCyclical=true ;
                fprintf(' Active-set is cyclical(#cycles=%i). \n',iCountCyclicActiveSet)
            end
            
            % I now have a dilemma, since the set has become cyclical it is
            % clear that if I deactivate any new nodes the thickness at
            % those nodes will become too small in the next active-set
            % interation. A solution is simply not to deactivate and to add
            % the deactivated nodes to the active set.
            
            if isActiveSetCyclical
                BCs1.hPosNode=[BCs1.hPosNode;Released] ; BCs1.hPosValue=BCs1.hPosNode*0+CtrlVar.ThickMin;
            end
            
            
        end
        
        
        
        if ~isActiveSetCyclical
            iCountCyclicActiveSet=0;
        end
        %%
        
        
        
    end   % active set loop
    
    
    
    
    
    %RunInfo.Iterations=nlIt(1); % here I use the number of NR iterations in the first update as a measure of
    % the nonlinearity of the problem
    
    % To do: I could consider relaxing the convergence criterias while
    % within the active set loop, and only fully converge once the active
    % set have been found, of course this assumes that as the final
    % convergence is reached, the active set no longer changes.
    
    if it > CtrlVar.ThicknessConstraintsItMax
        RunInfo.Forward.ActiveSetConverged=0;
    end
    
    
    if  it> CtrlVar.ThicknessConstraintsItMax
        fprintf(CtrlVar.fidlog,' Warning: In enforcing thickness constraints and finding a critical point, the loop was exited due to maximum number of iterations (%i) being reached. \n',CtrlVar.ThicknessConstraintsItMax);
    else
        if numel(BCs1.hPosNode)>0
            fprintf(CtrlVar.fidlog,'----- Active-set iteration converged after %-i iterations, constraining %-i thicknesses  \n \n ',it-1,numel(BCs1.hPosNode));
        end
    end
    
    
    if any(F1.h<CtrlVar.ThickMin)
        
        % Due to numerical errors it can sometimes happen that
        % on return even nodes in the active set
        % violate the pos thickness constraint.
        % If violated by less than 1e-10, I reset those to CtrlVar.ThickMin
        
        I=find(F1.h<CtrlVar.ThickMin) ;
        if max(F1.h(I))>(CtrlVar.ThickMin-1e-10)
            F1.h(I)=CtrlVar.ThickMin;
            if CtrlVar.ThicknessConstraintsInfoLevel>=1
                fprintf('Active-set: Resetting to limit. \n')
            end
        end
    end
    
    if any(F1.h<CtrlVar.ThickMin)
        save TestSave ;
        warning('some h1 <ThickMin on return from FIuvh2D. min(h1)=%-g',min(F1.h)) ;
        I=find(F1.h<CtrlVar.ThickMin) ;
        fprintf('Nodes with thickness<ThickMin: ') ; fprintf('%i ',I)  ; fprintf('\n')
        fprintf('                    thickness: ') ; fprintf('%g ',F1.h(I))  ; fprintf('\n')
        F1.h(F1.h<CtrlVar.ThickMin)=CtrlVar.ThickMin;
        %fprintf(CtrlVar.fidlog,' Found %-i thickness values less than %-g. Min thickness is %-g.',numel(indh0),CtrlVar.ThickMin,min(h));
        fprintf(CtrlVar.fidlog,' Setting h1(h1<%-g)=%-g \n ',CtrlVar.ThickMin,CtrlVar.ThickMin) ;
    end
end


end