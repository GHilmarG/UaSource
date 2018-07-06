function [ub1,vb1,ud1,vd1,h1,uvLambda,hLambdah,nlInfo,CtrlVar]=...
    FIuvh2D(CtrlVar,MUA,BCs,dt,S,B,ub0,vb0,ud0,vd0,h0,ub1,vb1,ud1,vd1,h1,as0,ab0,as1,ab1,...
    dubdt,dvbdt,duddt,dvddt,uvLambda,hLambdah,AGlen,C,n,m,alpha,rho,rhow,g)
    
    %        0  : values at t
    %        1  : on input (explicit) guess for values at t+dt, on output : converged values at t+dt
    %
    %

        
    if ~CtrlVar.ThicknessConstraints
        
        [ub1,vb1,ud1,vd1,h1,uvLambda,hLambdah,nlInfo]=uvh2D(CtrlVar,MUA,BCs,dt,h0,S,B,ub0,vb0,ud0,vd0,ub1,vb1,ud1,vd1,h1,as0,ab0,as1,ab1,dubdt,dvbdt,uvLambda,hLambdah,...
        AGlen,C,n,m,alpha,rho,rhow,g);
        
 
        
        CtrlVar.NumberOfActiveThicknedssConstraints=0;
        
    else
        %    NodesFixed: holdes the nodal numbers nodes in the active set
        %      ihactive: number nodes in active set
        %
        % -If min thickness is significanlty greater than CtrlVar.ThickMin, eliminate all thickness constraints
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
        ActiveSetReset=0;
        VariablesReset=0;
        hfixednodeposNew=[]; 
        if CtrlVar.ThicknessConstraintsInfoLevel>=10 ;
            fprintf(CtrlVar.fidlog,'  Enforcing min thickness of %-g using active-set method \n',CtrlVar.ThickMin);
        end
        
        if min(h1) > 1.1*CtrlVar.ThickMin;
            if CtrlVar.ThicknessConstraintsInfoLevel>=10 ;
                fprintf(CtrlVar.fidlog,' Eliminating any possible previous thickness constraints as min(h1)=%-g>1.1*CtrlVar.ThickMin=%-g \n',min(h0),CtrlVar.ThickMin);
            end
            Lhpos=[]; Lhrhspos=[] ; lambdahpos=[];
        end
        
        % possibly all thickness constraints were eliminated in AdapMesh when deactivating/activating elements
        % so I check if there are no thickness constrains but h0 is at the min thick.
        % if so then I introduce initial active set based on h0
        if isempty(Lhpos)
             Lhrhspos=[] ; lambdahpos=[];
            I=h1<=CtrlVar.ThickMin;
            Active=find(I);
            iActiveConstraints=numel(Active);
            if iActiveConstraints>0
                if CtrlVar.ThicknessConstraintsInfoLevel>=1 ;
                    fprintf(CtrlVar.fidlog,' Introducing %-i new initial active constraints based on h0 and h1  \n', iActiveConstraints);
                end
                [Lhpos,Lhrhspos]=createLh(MUA.Nnodes,Active,Active*0+CtrlVar.ThickMin,[],[]);
                lambdahpos=Lhrhspos*0; %  zeros(numel(Lhrhspos),1);
            end
            ub1(I)=0 ; vb1(I)=0; h1(I)=CtrlVar.ThickMin;
        end
        
        
        %         if any(h1<=CtrlVar.ThickMin)
        %             % not sure I fully understand why, but when new constraints must be introduced it
        %             % is for some reason better to reset initial estimates of the lagrange parameters
        %             %u1=u1*0 ; v1=v1*0 ;  % usually not needed, so try initially not to
        %             fprintf(CtrlVar.fidlog,' Resetting estimates for h1 and lambdahpos \n');
        %             h1=h0; lambdah=lambdah*0 ; lambdahpos=lambdahpos*0;
        %             u1=u0 ; v1=v0;
        %         end
        %
        
        isActiveSetModified=1;
        it=0;
        ItMax=CtrlVar.ThicknessConstraintsItMax  ;
        nlIt=zeros(ItMax+1,1);  %
        ihfixed=numel(hLambdah) ;          % these are the number of BCs applied to h by user, ie not including automatically added thickness constraints
        ihactive=numel(lambdahpos) ;      % number of active thickness constraints
        
        % iuvfixed=numel(lambdauv) ;        % these are the number of BCs applied to uv by user, ie not including automatically added thickness constraints
        % lambdauvpos=[];
        
        
        
        while isActiveSetModified==1  && it <= ItMax
            
            it=it+1;
            isActiveSetModified=0;
            
            if CtrlVar.ThicknessConstraintsInfoLevel>=1 ;
                if ihactive>0 || CtrlVar.ThicknessConstraintsInfoLevel>=10
                    fprintf(CtrlVar.fidlog,' Number of active thickness constraints is %-i \n',ihactive);
                end
            end
            
            %  if ~isempty(lambdahpos)   % I don't think I need this
            M=size(Lhpos,1); NodesFixed=zeros(M,1);
            for I=1:M
                NodesFixed(I)=find(Lhpos(I,:)==1);
            end
            ihactive=numel(NodesFixed) ;
                                    
            % identify elements where all thick nodal values are <=ThickMin
            % and set u1=0 ; v1=0; h1=ThickMin ;
            
            h1(h0<=CtrlVar.ThickMin)=CtrlVar.ThickMin;
            II=h1<=CtrlVar.ThickMin;
            h1(II)=CtrlVar.ThickMin;  ub1(II)=ub0(II) ; vb1(II)=vb0(II) ;  % modify initial guess for h1, important for convergence
            
            
            %             Setting mass-balance to zero is generally NOT a good idea.
            %             III=h1<=CtrlVar.ThickMin &   h0<=CtrlVar.ThickMin & (as1+ab1) <0 ;
            %             as0(III)=0 ; ab0(III)=0 ;  as1(III)=0 ; ab1(III)=0 ;
            %             if CtrlVar.ThicknessConstraintsInfoLevel>=10 ;
            %                 fprintf(CtrlVar.fidlog,' Mass balance set to zero at %-i nodes where balance is negative and thickness at ThickMin \n',numel(find(III)));
            %             end
            %
            
            uvhIt=1;
            
            while uvhIt<3
                % I have a possible inner solution here because
                % if uvh2D does not converge fully, it is usually best to just update the
                % active set on the partially converged solution
             
                Lh2=[Lh;Lhpos] ; Lhrhs2=[Lhrhs;Lhrhspos] ; lambdah2=[hLambdah;lambdahpos];
                
                
                [ub1,vb1,ud1,vd1,h1,uvLambda,lambdah2,nlInfo]=...
                    uvh2D(CtrlVar,MUA,dt,h0,S,B,ub0,vb0,ud0,vd0,ub1,vb1,ud1,vd1,h1,as0,ab0,as1,ab1,dubdt,dvbdt,Luv,Luvrhs,uvLambda,Lh2,Lhrhs2,lambdah2,...
                    AGlen,C,n,m,alpha,rho,rhow,g);
                
                
                if isempty(lambdah2)
                    hLambdah=[]; lambdahpos=[];
                else
                    lambdahpos=lambdah2(ihfixed+1:end); hLambdah=lambdah2(1:ihfixed);
                end
                
                
                if nlInfo.converged==1
                    %ActiveSetResetted=0;
                    VariablesReset=0;
                    break
                end

                
                if ~ActiveSetReset
                    ActiveSetReset=1;
                    if CtrlVar.ThicknessConstraintsInfoLevel>=1 ;
                        fprintf(CtrlVar.fidlog,' uvh2D did not converge in first active-set iteration. Eliminate active-set and try again \n');
                    end
                    Lhpos=[]; Lhrhspos=[] ; lambdahpos=[];  % eliminating all active constraints
                    ub1=ub0 ; vb1=vb0; h1=h0; 
                    NodesFixed=[];
                    ihactive=numel(NodesFixed) ;
                    
                elseif ~VariablesReset
                    VariablesReset=1;
                    if CtrlVar.ThicknessConstraintsInfoLevel>=1 ;
                        fprintf(CtrlVar.fidlog,' uvh2D did not converge. Resetting (u1,v1) to zero and setting h1=h0 \n');
                    end
                    ub1=ub0*0 ; vb1=vb0*0; h1=h0; hLambdah=hLambdah*0 ; lambdahpos=lambdahpos*0;
                end
                
                
                
                isActiveSetModified=1;
                uvhIt=uvhIt+1;
            end
            
            
            nlIt(it)=nlInfo.Iterations;  % if the active set is repeatadly updated, keep track of the number of nonlin iteration in each update
            
            
            %            lambdauvpos=lambdauv2(iuvfixed+1:end); lambdauv=lambdahuv2(1:iuvfixed);
            
            
            
            Info.NodesThickFixed=NodesFixed;
            Info.LambdaThickFixed=lambdahpos ;
            if CtrlVar.ThicknessConstraintsInfoLevel>=10 ;
                [~,I]=sort(lambdahpos);  % print out fixed nodes in the order of increasing lambda values
                fprintf(CtrlVar.fidlog,'            Nodes fixed: ')   ;
                fprintf(CtrlVar.fidlog,' \t %7i \t %7i \t %7i \t %7i \t %7i \t %7i \t %7i \t %7i \t %7i \t %7i \n \t \t \t \t \t \t',NodesFixed(I));
                fprintf(CtrlVar.fidlog,'\n   Lagrange multipliers: ') ;
                fprintf(CtrlVar.fidlog,' \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \t %+9.0f \n \t \t \t \t \t \t',lambdahpos(I)) ;
                fprintf(CtrlVar.fidlog,'\n');
            end
            
            
            
            
            
            % if the active-set method is selected, update active set
            % The active set is created/modified and the problem solved again if the active set has changed
            
            
            
            % Do I need to inactivate some thickness constraints?
            % if any of the lambdahpos are positive, then these constraints must be inactivated
            if ihactive>0   % are there any thickness constraints? If so see if some should be inactivated
                
                % sometimes constraints are being activated and inactivated over and over again. A possible remedy is not to inactivate constraints
                % immediately and to introduce a 1% threshold value
                
                lambdahposThreshold=CtrlVar.ThicknessConstraintsLambdaPosThreshold;
                if it>4 && CtrlVar.ThicknessConstraintsLambdaPosThreshold==0
                    
                    lambdahposThreshold=-mean(lambdahpos(lambdahpos<0))/100;
                    if isnan(lambdahposThreshold) ; lambdahposThreshold=0 ; end
                    if CtrlVar.ThicknessConstraintsInfoLevel>=1 ;
                        fprintf(CtrlVar.fidlog,' Introducing a min threshold of  %-g for  inactivating thickness constraints. \n',lambdahposThreshold);
                    end
                    
                end
                
                
                I=lambdahpos>lambdahposThreshold  ;  % if any of the Lagrange multipliers `lambdahpos' are positive, then these should be inactivated
                
                NewInActiveConstraints=find(I);
                
                iNewInActiveConstraints=numel(NewInActiveConstraints);
                if iNewInActiveConstraints>0   % have any become inactive?
                    Lhpos(I,:)=[]; Lhrhspos(I)=[] ; lambdahpos(I)=[];
                    isActiveSetModified=1;
                end
                if numel(Lhpos)==0  % if no constraints are left, reset to empty
                    Lhpos=[] ; Lhrhspos=[] ; lambdahpos=[];
                end
            else
                NewInActiveConstraints=[];
                iNewInActiveConstraints=numel(NewInActiveConstraints);
            end
            NodesReleased=NodesFixed(NewInActiveConstraints);
            
            NodesFixedOld=NodesFixed;
            M=size(Lhpos,1); NodesFixed=zeros(M,1);
            for I=1:M
                NodesFixed(I)=find(Lhpos(I,:)==1);
            end
            
            
            % Do I need to activate some new thickness constraints?
            % if any of h are less then MeshSizeMin then new constraints must be activated
            I=h1<=CtrlVar.ThickMin; % if thickness is less than ThickMin then further new thickness constraints must be introduced
            
            NodesWithTooSmallThick=find(I);
            NewActive=setdiff(NodesWithTooSmallThick,NodesFixed);
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
            hfixednodeposNew=NewActive ;
            if iNewActiveConstraints>0
                if CtrlVar.ThicknessConstraintsInfoLevel>=1 ;
                    fprintf(CtrlVar.fidlog,' %i new active constraints \n',iNewActiveConstraints);
                end
                isActiveSetModified=1;
                hfixednodeposNew=NewActive ; hfixedvalueposNew=hfixednodeposNew*0+CtrlVar.ThickMin;
                [LhposNew,LhrhsposNew]=createLh(MUA.Nnodes,hfixednodeposNew,hfixedvalueposNew,[],[]);
                lambdahposNew=LhrhsposNew*0 ; %   zeros(numel(LhrhsposNew),1);
                Lhpos=[Lhpos;LhposNew] ; Lhrhspos=[Lhrhspos;LhrhsposNew] ; lambdahpos=[lambdahpos;lambdahposNew];
                
            end
            
            M=size(Lhpos,1); NodesFixed=zeros(M,1);
            for I=1:M
                NodesFixed(I)=find(Lhpos(I,:)==1);
            end
            
            %setdiff(NodesFixed,NodesFixedOld)
            
            ihactive=numel(NodesFixed) ;
            ihfixed=numel(hLambdah) ;
            
            
            % modify initial guess for h1, important for convergence
            %h1(NewActive)=ThickMin;
            
            if CtrlVar.ThicknessConstraintsInfoLevel>=1 ;
                if iNewInActiveConstraints> 0 || iNewActiveConstraints> 0
                    fprintf(CtrlVar.fidlog,' Updating thickness constraints: inactivated: %-i,  activated: %-i, total number of thickness constrains: %-i \n',...
                        iNewInActiveConstraints,iNewActiveConstraints,ihactive);
                    fprintf(CtrlVar.fidlog,'  Nodes inactivated: ')   ;
                    fprintf(CtrlVar.fidlog,' \t %7i \t %7i \t %7i \t %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \n \t \t \t \t \t',NodesReleased);
                    
                    fprintf(CtrlVar.fidlog,'\n    Nodes activated: ')   ;
                    fprintf(CtrlVar.fidlog,' \t %7i \t %7i \t %7i \t %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \t  %7i \n \t \t \t \t \t',hfixednodeposNew);
                    fprintf(CtrlVar.fidlog,'\n ')   ;
                    
                    
                end
                
                if isActiveSetModified==1 && it <= ItMax
                    fprintf(CtrlVar.fidlog,' Active set modifed. System is solved again using the new active set. ActiveSet Iteration Nr. %-i \n',it);
                end
            end
        end
        
        %nlInfo.Iterations=nlIt(1); % here I use the number of NR iterations in the first update as a measure of
        % the nonlinearity of the problem
        if CtrlVar.ThicknessConstraintsInfoLevel>=1 ;
            if isActiveSetModified==1  &&  it> ItMax
                fprintf(CtrlVar.fidlog,' Warning: In enforcing thickness constraints and finding a critical point, the loop was exited due to maximum number of iterations being reached. \n');
            end
            
            if isActiveSetModified==0 && ihactive>0
                fprintf(CtrlVar.fidlog,'===== ActiveSet iteration converged after %-i iterations, constraining %-i thicknesses  \n \n ',it-1,ihactive);
            end
        end
        
        CtrlVar.NumberOfActiveThicknessConstraints=ihactive;
    end
   
    Info.nlInfo=nlInfo;
    
end