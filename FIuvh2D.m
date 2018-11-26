function [UserVar,RunInfo,F1,l1,BCs1,dt]=FIuvh2D(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1)


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
    ItMax=CtrlVar.ThicknessConstraintsItMax  ;
    nlIt=zeros(ItMax+1,1)*2+NaN;  %
    Released=[] ; Activated=[];
    
    while isActiveSetModified==1
        
        it=it+1;
        isActiveSetModified=0;
        
        if CtrlVar.ThicknessConstraintsInfoLevel>=1 
            if numel(BCs1.hPosNode)>0 || CtrlVar.ThicknessConstraintsInfoLevel>=10
                fprintf(CtrlVar.fidlog,' Number of active thickness constraints is %-i \n',numel(BCs1.hPosNode));
            end
        end
        
        II= (F0.h<=CtrlVar.ThickMin) | (F1.h<=CtrlVar.ThickMin) ;
        F1.h(II)=CtrlVar.ThickMin;  F1.ub(II)=F0.ub(II) ; F1.vb(II)=F0.vb(II) ;  % modify initial guess for h1, important for convergence
        
        F1.h(BCs1.hPosNode)=CtrlVar.ThickMin;
        
        
        uvhIt=0; ReduceTimeStep=0;% needed in inner loop
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
            lambdahpos=hLambda(numel(BCs1.hFixedNode)+numel(BCs1.hTiedNodeA)+1:end) ;%  if I always put the hPos constraints at end of all other h constraints
            % then this will work
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
                    fprintf(CtrlVar.fidlog,' uvh2D did not converge. Resetting u1, v1 and h1 to values at start of time step. \n');
                end
                
                
                if CtrlVar.WriteRunInfoFile
                    fprintf(RunInfo.File.fid,' uvh2D did not converge. Resetting u1, v1 and h1 to values at start of time step. \n');
                end
                
                
                F1.ub=F0.ub ; F1.vb=F0.vb; F1.ud=F0.ud ; F1.vd=F0.vd ; F1.h=F0.h;
                F1.h(BCs1.hPosNode)=CtrlVar.ThickMin;  % consider adding
               
                l1.ubvb=l1.ubvb*0 ; l1.udvd=l1.udvd*0; l1.h=l1.h*0;
               % BCs1.hPosNode=[] ;  BCs1.hPosValue=[]; isActiveSetModified=1;
%             elseif ~ActiveSetReset
%                 
%                 ActiveSetReset=1;
%                 if CtrlVar.ThicknessConstraintsInfoLevel>=1
%                     fprintf(CtrlVar.fidlog,' uvh2D did not converge in first active-set iteration. Resetting variables and eliminating active-set. \n');
%                 end
%                 F1.ub=F0.ub ; F1.vb=F0.vb; F1.ud=F0.ud ; F1.vd=F0.vd ; F1.h=F0.h; 
%         
%                 l1.ubvb=l1.ubvb*0 ; l1.udvd=l1.udvd*0; l1.h=l1.h*0;
%                 BCs1.hPosNode=[] ;  BCs1.hPosValue=[];
%                 isActiveSetModified=1;
%
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
            
            
            
            
           
        end
        
        
        
        % keep track of the number of non-lin iteration in each update
        
        
        if CtrlVar.ThicknessConstraintsInfoLevel>=10 
            [~,I]=sort(lambdahpos);  % print out fixed nodes in the order of increasing lambda values
            fprintf(CtrlVar.fidlog,'            Nodes fixed: ')   ;
            fprintf(CtrlVar.fidlog,' \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \n \t \t \t \t \t \t',BCs1.hPosNode(I));
            fprintf(CtrlVar.fidlog,'\n   Lagrange multipliers: ') ;
            fprintf(CtrlVar.fidlog,' \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \n \t \t \t \t \t \t',hLambda(I)) ;
            fprintf(CtrlVar.fidlog,'\n');
        end
        
        
        if  it > ItMax
            break   % I must exit before the active set is (potentally) modified further, or the Lagrange variable
            % calculated in the call to uvh2D may not be consistent with the active set.
        end
        
        
        
        % Now I've solved for h1 and if needed a new active set must be defined
        %
        % The new active set contains all nodes where h1 less than hmin that were not in the previous active set
        % Those of the nodes in the previous active set with positve slack values
        % Nodes in the previous set with negative slack values must be taken out of the set
        
        
        % if the active-set method is selected, update active set
        % The active set is created/modified and the problem solved again if the active set has changed
        
        
        % Do I need to inactivate some thickness constraints?
        % if any of the lambdahpos are positive, then these constraints must be inactivated
        
        if numel(BCs1.hPosNode)>0   % are there any thickness constraints? If so see if some should be inactivated
            
            % sometimes constraints are being activated and inactivated over and over again. A possible remedy is not to inactivate constraints
            % immediately and to introduce a 1% threshold value
            
            lambdahposThreshold=CtrlVar.ThicknessConstraintsLambdaPosThreshold;
%             
%             if it>4 && CtrlVar.ThicknessConstraintsLambdaPosThreshold==0
%                 
%                 lambdahposThreshold=-mean(lambdahpos(lambdahpos<0))/100;
%                 if isnan(lambdahposThreshold) ; lambdahposThreshold=0 ; end
%                 if CtrlVar.ThicknessConstraintsInfoLevel>=1 ;
%                     fprintf(CtrlVar.fidlog,' Introducing a min threshold of  %-g for  inactivating thickness constraints. \n',lambdahposThreshold);
%                 end
%             end
            
            
            I=lambdahpos>lambdahposThreshold  ;  % if any of the Lagrange multipliers `lambdahpos' are positive, then these should be inactivated
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
                fprintf(CtrlVar.fidlog,' Number of new active thickness constraints %-i larger then max number or newly added constraints %-i \n ',...
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
        
        LastReleased=Released; LastActivated=Activated;
        Released=setdiff(LastActiveSet,BCs1.hPosNode)   ;% nodes in last active set that are no longer in the new one
        Activated=setdiff(BCs1.hPosNode,LastActiveSet)  ;% nodes in new active set that were not in the previous one
        
        if ~isempty(LastReleased)
            if isempty(setxor(LastReleased,Activated))
                fprintf(' Previous in-activated node set equal to the currently activated one. \n')
            end
        end
        
        if ~isempty(LastActivated)
            if isempty(setxor(LastActivated,Released))
                fprintf(' Previous activated node set equal to the currently in-activated one. \n')
            end
        end
        
        if CtrlVar.ThicknessConstraintsInfoLevel>=1 
            if iNewInActiveConstraints> 0 || iNewActiveConstraints> 0
                fprintf(CtrlVar.fidlog,' Updating thickness constraints: inactivated: %-i,  activated: %-i, total number of thickness constrains: %-i \n',...
                    iNewInActiveConstraints,iNewActiveConstraints,numel(BCs1.hPosNode));
                fprintf(CtrlVar.fidlog,'  Nodes inactivated: ')   ;
                fprintf(CtrlVar.fidlog,' \t %7i \t %7i \t %7i \t %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \n \t \t \t \t \t',NodesReleased);
                
                fprintf(CtrlVar.fidlog,'\n    Nodes activated: ')   ;
                fprintf(CtrlVar.fidlog,' \t %7i \t %7i \t %7i \t %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \n \t \t \t \t \t',NewActive);
                fprintf(CtrlVar.fidlog,'\n ')   ;
            else
                fprintf(CtrlVar.fidlog,'No thickness constraints activated or deactivated. \n')   ;
            end
            
            if isActiveSetModified==1 && it <= ItMax
                fprintf(CtrlVar.fidlog,' Active set modified. System is solved again using the new active set. ActiveSet Iteration Nr. %-i \n',it);
            end
        end
    end
    
    
    
    
    
    %RunInfo.Iterations=nlIt(1); % here I use the number of NR iterations in the first update as a measure of
    % the nonlinearity of the problem
    
    % To do: I could consider relaxing the convergence criterias while
    % within the active set loop, and only fully converge once the active
    % set have been found, of course this assumes that as the final
    % convergence is reached, the active set no longer changes.
    
    if it > ItMax
        RunInfo.Forward.ActiveSetConverged=0;
    end
    
    
    if  it> ItMax
        fprintf(CtrlVar.fidlog,' Warning: In enforcing thickness constraints and finding a critical point, the loop was exited due to maximum number of iterations (%i) being reached. \n',ItMax);
    else
        if numel(BCs1.hPosNode)>0
            fprintf(CtrlVar.fidlog,'----- ActiveSet iteration converged after %-i iterations, constraining %-i thicknesses  \n \n ',it-1,numel(BCs1.hPosNode));
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
            if CtrlVar.ThicknessConstraintsInfoLevel>=10 
                fprintf('Resetting to limit\n')
            end
        end
    end
    
    if any(F1.h<CtrlVar.ThickMin) 
        save TestSave ;
        warning('some h1 <ThickMin on return from FIuvh2D. min(h1)=%-g',max(F1.h)) ;
        I=find(F1.h<CtrlVar.ThickMin) ;
        fprintf('Nodes with thickness<ThickMin: ') ; fprintf('%i ',I)  ; fprintf('\n')
        fprintf('                    thickness: ') ; fprintf('%g ',F1.h(I))  ; fprintf('\n')
        F1.h(F1.h<CtrlVar.ThickMin)=CtrlVar.ThickMin;
        %fprintf(CtrlVar.fidlog,' Found %-i thickness values less than %-g. Min thickness is %-g.',numel(indh0),CtrlVar.ThickMin,min(h));
        fprintf(CtrlVar.fidlog,' Setting h1(h1<%-g)=%-g \n ',CtrlVar.ThickMin,CtrlVar.ThickMin) ;
    end
end


end