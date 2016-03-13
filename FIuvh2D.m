function [ub1,vb1,ud1,vd1,h1,uvLambda,hLambda,RunInfo,CtrlVar,BCs,dt]=...
    FIuvh2D(CtrlVar,MUA,BCs,dt,S,B,ub0,vb0,ud0,vd0,h0,ub1,vb1,ud1,vd1,h1,as0,ab0,as1,ab1,...
    dubdt,dvbdt,duddt,dvddt,uvLambda,hLambda,AGlen,C,n,m,alpha,rho,rhow,g)

%        0  : values at t
%        1  : on input (explicit) guess for values at t+dt, on output : converged values at t+dt
%
%

RunInfo.ActiveSetConverged=1;

if ~CtrlVar.ThicknessConstraints
    
    [ub1,vb1,ud1,vd1,h1,uvLambda,hLambda,RunInfo]=uvh2D(CtrlVar,MUA,BCs,dt,h0,S,B,ub0,vb0,ud0,vd0,ub1,vb1,ud1,vd1,h1,as0,ab0,as1,ab1,dubdt,dvbdt,uvLambda,hLambda,...
        AGlen,C,n,m,alpha,rho,rhow,g);
    
    
    
    CtrlVar.NumberOfActiveThicknedssConstraints=0;
    
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
    
    
    
    
    %hfixednodeposNew=[];
    if CtrlVar.ThicknessConstraintsInfoLevel>=10 ;
        fprintf(CtrlVar.fidlog,'  Enforcing min thickness of %-g using active-set method \n',CtrlVar.ThickMin);
    end
    
    if min(h1) > 1.1*CtrlVar.ThickMin;
        if CtrlVar.ThicknessConstraintsInfoLevel>=10 ;
            fprintf(CtrlVar.fidlog,' Eliminating any possible previous thickness constraints as min(h1)=%-g>1.1*CtrlVar.ThickMin=%-g \n',min(h0),CtrlVar.ThickMin);
        end
        BCs.hPosNode=[] ; BCs.hPosValue=[];
    end
    
    % possibly all thickness constraints were eliminated in AdapMesh when deactivating/activating elements
    % so I check if there are no thickness constrains but h0 is at the min thick.
    % if so then I introduce an initial active set based on h0
    %!if isempty(Lhpos)
    if isempty(BCs.hPosNode)
        Active=find(h0<=CtrlVar.ThickMin);
        BCs.hPosNode=Active ; BCs.hPosValue=BCs.hPosNode*0+CtrlVar.ThickMin;
        ub1(BCs.hPosNode)=0 ; vb1(BCs.hPosNode)=0; h1(BCs.hPosNode)=CtrlVar.ThickMin; % now set estimates of velocities at hPosNode to zero
        if numel(BCs.hPosNode)>0
            if CtrlVar.ThicknessConstraintsInfoLevel>=1 ;
                fprintf(CtrlVar.fidlog,' Introducing %-i new initial active constraints based on h0  \n', numel(BCs.hPosNode));
            end
        end
    end
    
    isActiveSetModified=1;
    it=0;
    ItMax=CtrlVar.ThicknessConstraintsItMax  ;
    nlIt=zeros(ItMax+1,1);  %
    Released=[] ; Activated=[];
    
    while isActiveSetModified==1
        
        it=it+1;
        isActiveSetModified=0;
        
        if CtrlVar.ThicknessConstraintsInfoLevel>=1 ;
            if numel(BCs.hPosNode)>0 || CtrlVar.ThicknessConstraintsInfoLevel>=10
                fprintf(CtrlVar.fidlog,' Number of active thickness constraints is %-i \n',numel(BCs.hPosNode));
            end
        end
        
        II= (h0<=CtrlVar.ThickMin) | (h1<=CtrlVar.ThickMin) ;
        h1(II)=CtrlVar.ThickMin;  ub1(II)=ub0(II) ; vb1(II)=vb0(II) ;  % modify initial guess for h1, important for convergence
        
        h1(BCs.hPosNode)=CtrlVar.ThickMin;
        
        
        uvhIt=1; ActiveSetReset=0; VariablesReset=0;  ReduceTimeStep=0;% needed in inner loop
        while uvhIt<3  % I have a possible inner iteration here because
            % if uvh2D does not converge fully, it is usually best to just update the
            % active set on the partially converged solution
            
            %if any(h0<CtrlVar.ThickMin) ;
            %    save TestSave ; warning('h0 negative on input to uvh2D. min(h0)=%g ') ;
            %    keyboard
            %end
            
            [ub1,vb1,ud1,vd1,h1,uvLambda,hLambda,RunInfo]=uvh2D(CtrlVar,MUA,BCs,dt,h0,S,B,ub0,vb0,ud0,vd0,ub1,vb1,ud1,vd1,h1,as0,ab0,as1,ab1,dubdt,dvbdt,uvLambda,hLambda,...
                AGlen,C,n,m,alpha,rho,rhow,g);
            
            
            
            % keep a copy of the old active set
            LastActiveSet=BCs.hPosNode;
            
            %must now find the lambda values corresponding to nodes that where constrainted to pos thickness
            
            lambdahpos=hLambda(numel(BCs.hFixedNode)+1:end);  % if I always put the hPos constraints at end of all other h constraints
            % then this will work
            if numel(lambdahpos) ~= numel(BCs.hPosNode)
                save TestSave ; error(' # of elements in lambdahpos must equal # of elements in BCs.hPosNode')
            end
            
            if RunInfo.converged==1
                break
            end
            
            if ~ReduceTimeStep
                ReduceTimeStep=1;
                fprintf(CtrlVar.fidlog,' Warning : Reducing time step from %-g to %-g \n',dt,dt/10);
                dt=dt/10; CtrlVar.dt=dt;
                fprintf(CtrlVar.fidlog,'Also resetting u1, v1, h1 to ub0, vb0 and h0, and setting estimates for Lagrange parameters to zero. \n');
                ub1=ub0*0 ; vb1=vb0*0 ;  ud1=ud0*0 ; vd1=vd0*0 ; h1=h0;
                uvLambda=uvLambda*0; hLambda=hLambda*0;
            elseif ~ActiveSetReset
                ActiveSetReset=1;
                if CtrlVar.ThicknessConstraintsInfoLevel>=1 ;
                    fprintf(CtrlVar.fidlog,' uvh2D did not converge in first active-set iteration. Eliminate active-set and try again \n');
                end
                ub1=ub0 ; vb1=vb0; h1=h0; uvLambda=uvLambda*0; hLambda=hLambda*0;
                BCs.hPosNode=[] ;  BCs.hPosValue=[] ;
                
            elseif ~VariablesReset
                VariablesReset=1;
                if CtrlVar.ThicknessConstraintsInfoLevel>=1 ;
                    fprintf(CtrlVar.fidlog,' uvh2D did not converge. Resetting (u1,v1) to zero and setting h1=h0 \n');
                end
                ub1=ub0*0 ; vb1=vb0*0; h1=h0; uvLambda=uvLambda*0; hLambda=hLambda*0;
            end
            
            isActiveSetModified=1;
            uvhIt=uvhIt+1;
        end
        
        
        nlIt(it)=RunInfo.Iterations;  % if the active set is repeatedly updated,
        % keep track of the number of non-lin iteration in each update
        
        
        if CtrlVar.ThicknessConstraintsInfoLevel>=10 ;
            [~,I]=sort(lambdahpos);  % print out fixed nodes in the order of increasing lambda values
            fprintf(CtrlVar.fidlog,'            Nodes fixed: ')   ;
            fprintf(CtrlVar.fidlog,' \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \t %9i \n \t \t \t \t \t \t',BCs.hPosNode(I));
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
        
        if numel(BCs.hPosNode)>0   % are there any thickness constraints? If so see if some should be inactivated
            
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
                BCs.hPosNode(I)=[]; BCs.hPosValue(I)=[];
                isActiveSetModified=1;
            end
            
        else
            NewInActiveConstraints=[];
            iNewInActiveConstraints=numel(NewInActiveConstraints);
        end
        
        NodesReleased=LastActiveSet(NewInActiveConstraints);
        
        % Do I need to activate some new thickness constraints?
        %I=h1<=CtrlVar.ThickMin; % if thickness is less than ThickMin then further new thickness constraints must be introduced
        I=h1<=(CtrlVar.ThickMin-100*eps); % if thickness is less than ThickMin then further new thickness constraints must be introduced
        
        NodesWithTooSmallThick=find(I);
        NewActive=setdiff(NodesWithTooSmallThick,BCs.hPosNode);  % exclude those already in the active set
        NewActive=setdiff(NewActive,NodesReleased);  % do not include those nodes at min thick that I now must release
        
        
        iNewActiveConstraints=numel(NewActive);
        
        if iNewActiveConstraints> CtrlVar.MaxNumberOfNewlyIntroducedActiveThicknessConstraints
            if CtrlVar.ThicknessConstraintsInfoLevel>=1 ;
                fprintf(CtrlVar.fidlog,' Number of new active thickness constraints %-i larger then max number or newly added constraints %-i \n ',...
                    iNewActiveConstraints,CtrlVar.MaxNumberOfNewlyIntroducedActiveThicknessConstraints);
                fprintf(CtrlVar.fidlog,' Only the smallest %-i thickness values are constrained \n',CtrlVar.MaxNumberOfNewlyIntroducedActiveThicknessConstraints);
            end
            [Temp,II]=sort(h1);
            NewActive=II(1:CtrlVar.MaxNumberOfNewlyIntroducedActiveThicknessConstraints);
            iNewActiveConstraints=CtrlVar.MaxNumberOfNewlyIntroducedActiveThicknessConstraints;
            
        end
        
        if iNewActiveConstraints>0
            if CtrlVar.ThicknessConstraintsInfoLevel>=1 ;
                fprintf(CtrlVar.fidlog,' %i new active constraints \n',iNewActiveConstraints);
            end
            isActiveSetModified=1;
            BCs.hPosNode=[BCs.hPosNode;NewActive] ; BCs.hPosValue=BCs.hPosNode*0+CtrlVar.ThickMin;
        end
        
        
        % modify initial guess for h1, important for convergence
        %h1(NewActive)=ThickMin;
        
        LastReleased=Released; LastActivated=Activated;
        Released=setdiff(LastActiveSet,BCs.hPosNode)   ;% nodes in last active set that are no longer in the new one
        Activated=setdiff(BCs.hPosNode,LastActiveSet)  ;% nodes in new active set that were not in the previous one
        
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
        
        if CtrlVar.ThicknessConstraintsInfoLevel>=1 ;
            if iNewInActiveConstraints> 0 || iNewActiveConstraints> 0
                fprintf(CtrlVar.fidlog,' Updating thickness constraints: inactivated: %-i,  activated: %-i, total number of thickness constrains: %-i \n',...
                    iNewInActiveConstraints,iNewActiveConstraints,numel(BCs.hPosNode));
                fprintf(CtrlVar.fidlog,'  Nodes inactivated: ')   ;
                fprintf(CtrlVar.fidlog,' \t %7i \t %7i \t %7i \t %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \n \t \t \t \t \t',NodesReleased);
                
                fprintf(CtrlVar.fidlog,'\n    Nodes activated: ')   ;
                fprintf(CtrlVar.fidlog,' \t %7i \t %7i \t %7i \t %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \n \t \t \t \t \t',NewActive);
                fprintf(CtrlVar.fidlog,'\n ')   ;
            else
                fprintf(CtrlVar.fidlog,'No nodes activated or deactivated. \n')   ;
            end
            
            if isActiveSetModified==1 && it <= ItMax
                fprintf(CtrlVar.fidlog,' Active set modifed. System is solved again using the new active set. ActiveSet Iteration Nr. %-i \n',it);
            end
        end
    end
    
    
    
    
    
    %RunInfo.Iterations=nlIt(1); % here I use the number of NR iterations in the first update as a measure of
    % the nonlinearity of the problem
    
    if it > ItMax
        RunInfo.ActiveSetConverged=0;
    end
    
    
    if  it> ItMax
        fprintf(CtrlVar.fidlog,' Warning: In enforcing thickness constraints and finding a critical point, the loop was exited due to maximum number of iterations (%i) being reached. \n',ItMax);
    else
        if numel(BCs.hPosNode)>0
            fprintf(CtrlVar.fidlog,'----- ActiveSet iteration converged after %-i iterations, constraining %-i thicknesses  \n \n ',it-1,numel(BCs.hPosNode));
        end
    end
    
    
    if any(h1<CtrlVar.ThickMin)
        
        % Due to numerical errors it can sometimes happen that
        % on return even nodes in the active set
        % violate the pos thickness constraint.
        % If violated by less than 1e-10, I reset those to CtrlVar.ThickMin
        
        I=find(h1<CtrlVar.ThickMin) ;
        if max(h1(I))>(CtrlVar.ThickMin-1e-10) ;
            h1(I)=CtrlVar.ThickMin;
            if CtrlVar.ThicknessConstraintsInfoLevel>=10 ;
                fprintf('Resetting to limit\n')
            end
        end
    end
    
    if any(h1<CtrlVar.ThickMin) ;
        save TestSave ;
        warning('some h1 <ThickMin on return from FIuvh2D. min(h1)=%-g',max(h1)) ;
        I=find(h1<CtrlVar.ThickMin) ;
        fprintf('Nodes with thickness<ThickMin: ') ; fprintf('%i ',I)  ; fprintf('\n')
        fprintf('                    thickness: ') ; fprintf('%g ',h1(I))  ; fprintf('\n')
        h1(h1<CtrlVar.ThickMin)=CtrlVar.ThickMin;
        %fprintf(CtrlVar.fidlog,' Found %-i thickness values less than %-g. Min thickness is %-g.',numel(indh0),CtrlVar.ThickMin,min(h));
        fprintf(CtrlVar.fidlog,' Setting h1(h1<%-g)=%-g \n ',CtrlVar.ThickMin,CtrlVar.ThickMin) ;
    end
end


end