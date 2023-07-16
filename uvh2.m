function [UserVar,RunInfo,F1,l1,BCs1,dt]=uvh2(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1)
    

    narginchk(9,9)
    nargoutchk(6,6)
    
    dt=CtrlVar.dt;

    RunInfo.Forward.ActiveSetConverged=1;
    RunInfo.Forward.uvhIterationsTotal=0;
    iActiveSetIteration=0;
    isActiveSetCyclical=NaN;

    if CtrlVar.LevelSetMethod &&  ~isnan(CtrlVar.LevelSetDownstreamAGlen) &&  ~isnan(CtrlVar.LevelSetDownstream_nGlen)
        
        if isempty(F0.LSFMask)  % If I have already solved the LSF equation, this will not be empty and does not need to be recalculated (ToDo)
            F0.LSFMask=CalcMeshMask(CtrlVar,MUA,F0.LSF,0);
        end

        F1.LSFMask=F0.LSFMask;
        if ~isnan(CtrlVar.LevelSetDownstreamAGlen)
            F0.AGlen(F0.LSFMask.NodesOut)=CtrlVar.LevelSetDownstreamAGlen;
            F1.AGlen(F1.LSFMask.NodesOut)=CtrlVar.LevelSetDownstreamAGlen;
            F0.n(F0.LSFMask.NodesOut)=CtrlVar.LevelSetDownstream_nGlen;
            F1.n(F1.LSFMask.NodesOut)=CtrlVar.LevelSetDownstream_nGlen;
        end

    end


    if ~CtrlVar.ThicknessConstraints


        [UserVar,RunInfo,F1,l1,BCs1]=uvh2D(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);

        if ~RunInfo.Forward.uvhConverged

            [UserVar,RunInfo,F1,F0,l0,l1,BCs1,dt]=uvh2NotConvergent(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1) ;

        end

        F1=UpdateFtimeDerivatives(UserVar,RunInfo,CtrlVar,MUA,F1,F0) ; % but currently F0 is not returned...


        if numel(RunInfo.Forward.uvhActiveSetIterations)<CtrlVar.CurrentRunStepNumber
            RunInfo.Forward.uvhActiveSetIterations=[RunInfo.Forward.uvhActiveSetIterations;RunInfo.Forward.uvhActiveSetIterations+NaN];
            RunInfo.Forward.uvhActiveSetCyclical=[RunInfo.Forward.uvhActiveSetCyclical;RunInfo.Forward.uvhActiveSetCyclical+NaN];
            RunInfo.Forward.uvhActiveSetConstraints=[RunInfo.Forward.uvhActiveSetConstraints;RunInfo.Forward.uvhActiveSetConstraints+NaN];
        end

        RunInfo.Forward.uvhActiveSetIterations(CtrlVar.CurrentRunStepNumber)=NaN ;
        RunInfo.Forward.uvhActiveSetCyclical(CtrlVar.CurrentRunStepNumber)=NaN;
        RunInfo.Forward.uvhActiveSetConstraints(CtrlVar.CurrentRunStepNumber)=NaN;






    else   %  Thickness constraints used

        [UserVar,RunInfo,F1,l1,BCs1,isActiveSetModified,Activated,Released]=ActiveSetInitialisation(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1) ;
        LastReleased=Released;
        LastActivated=Activated;


        VariablesReset=0;  % keep this out here, only set once, whereas ReduceTimeStep is reset in each active-set iteration.
        iCounter=1;
        
        % Released=[] ; Activated=[];
        
        
        
        while true   % active-set loop
            
    
            
            % iActiveSetIteration=iActiveSetIteration+1;
            
            
            fprintf(CtrlVar.fidlog,' ----------> Active-set Iteration #%-i <----------  \n',iActiveSetIteration);
            
            if CtrlVar.ThicknessConstraintsInfoLevel>=1
                if numel(BCs1.hPosNode)>0 || CtrlVar.ThicknessConstraintsInfoLevel>=10
                    fprintf(CtrlVar.fidlog,' Number of active thickness constraints is %-i \n',numel(BCs1.hPosNode));
                end
            end

            
            uvhIt=0; ReduceTimeStep=0;% needed in inner loop
            
            %% Solve uv for the current active set
            while uvhIt<=2  % I have a possible inner iteration here because
                
                uvhIt=uvhIt+1;
               
               
                

                if CtrlVar.DebugMode
                    warning('FIuvh2D:nonconvergent','Ahead of uvh2D call')
                    filename='Dump_Before_FIuvh2D.mat';
                    fprintf('Saving all data in %s \n',filename)
                    save(filename)
                end

                F1.h(BCs1.hPosNode)=CtrlVar.ThickMin;  % Make sure iterate is feasible, but this should be delt with by the BCs anyhow
                
                F1.ub(BCs1.hPosNode)=F0.ub(BCs1.hPosNode);  % might help with convergence
                F1.vb(BCs1.hPosNode)=F0.vb(BCs1.hPosNode);

                [UserVar,RunInfo,F1,l1,BCs1]=uvh2D(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);


                switch CtrlVar.FlowApproximation
                    case "SSTREAM"

                        nlIt(iCounter)=RunInfo.Forward.uvhIterations(CtrlVar.CurrentRunStepNumber);  iCounter=iCounter+1;

                        % so what is the best estimate for number of NR-uvh iterations over the active set iteration
                        % It is possible that if several active set iterations are done, the mean value will be
                        % so small that selecting a time step based on that mean value will cause convergence issues in the first active set
                        % iteration at a later stage.

                       % RunInfo.Forward.uvhIterations(CtrlVar.CurrentRunStepNumber)=mean(nlIt,'omitnan');
                        RunInfo.Forward.uvhIterations(CtrlVar.CurrentRunStepNumber)=max(nlIt,[],'omitnan');


                    case "SSHEET"
                        nlIt(iCounter)=RunInfo.Forward.hIterations(CtrlVar.CurrentRunStepNumber);  iCounter=iCounter+1;
                        RunInfo.Forward.hIterations(CtrlVar.CurrentRunStepNumber)=mean(nlIt,'omitnan');
                end


                if RunInfo.Forward.uvhConverged==1
                    break
                end


                warning('uvh:uvhSolutionNotConvergent','uvh2D did not converge')
                filename='Dump_uvh';
                fprintf('Saving all local data in %s \n',filename)
                try

                    save(filename)

                catch ME

                    disp('Error Message:')
                    disp(ME.message)
                    warning('uvh:CouldNotSaveFile','For some reason file could not be saved.')

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
                    dt=dt/2; CtrlVar.dt=dt;  F1.dt=CtrlVar.dt ;  F0.dt=CtrlVar.dt ; 
                    fprintf(CtrlVar.fidlog,' Warning : Reducing time step from %-g to %-g \n',dtOld,CtrlVar.dt);
                    
                    F0.h(F0.h<=CtrlVar.ThickMin)=CtrlVar.ThickMin ; F0.h(BCs1.hPosNode)=CtrlVar.ThickMin;
                    [UserVar,RunInfo,F0,l0]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs1,F0,l0);
           
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
                    dt=dt/2; CtrlVar.dt=dt; F1.dt=CtrlVar.dt ;  F0.dt=CtrlVar.dt ; 
                    
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



            % Update active set once uvh solution has been found.
            % Note: If this active-set iteration does not converge,
            %       the uvh solution and the active set are not consistence,
            %       i.e. the uvh solution was not obtained for the BCs in the active set.
            [UserVar,RunInfo,BCs1,lambdahpos,isActiveSetModified,isActiveSetCyclical,Activated,Released]=ActiveSetUpdate(UserVar,RunInfo,CtrlVar,MUA,F1,l1,BCs1,iActiveSetIteration,LastReleased,LastActivated);

            LastReleased=Released;
            LastActivated=Activated;


        

            if ~isActiveSetModified
                fprintf(' Leaving active-set loop because active set did not change in last active-set iteration. \n')
                break

            end

            if isActiveSetCyclical

                fprintf(" Leaving active-set pos. thickness loop because it has become cyclical\n");
                break

            end


            iActiveSetIteration=iActiveSetIteration+1;
            if  iActiveSetIteration > CtrlVar.ThicknessConstraintsItMax
                RunInfo.Forward.ActiveSetConverged=0;
                fprintf(' Leaving active-set pos. thickness loop because number of active-set iteration (%i) greater than maximum allowed (CtrlVar.ThicknessConstraintsItMax=%i). \n ',iActiveSetIteration,CtrlVar.ThicknessConstraintsItMax)
                break
            end



        end   % active set loop



        % To do: I could consider relaxing the convergence criterias while
        % within the active set loop, and only fully converge once the active
        % set have been found, of course this assumes that as the final
        % convergence is reached, the active set no longer changes.
        
     
        
        
        if  ~RunInfo.Forward.ActiveSetConverged
            fprintf(CtrlVar.fidlog,' Warning: In enforcing thickness constraints and finding a critical point, the loop was exited due to maximum number of iterations (%i) being reached. \n',CtrlVar.ThicknessConstraintsItMax);
        else
            if numel(BCs1.hPosNode)>0
                fprintf(CtrlVar.fidlog,'----- Active-set iteration converged after %-i iterations, constraining %-i thicknesses  \n \n ',iActiveSetIteration-1,numel(BCs1.hPosNode));
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

            warning('some h1 <ThickMin on return from FIuvh2D. min(h1)=%-g',min(F1.h)) ;
            I=find(F1.h<CtrlVar.ThickMin) ;
            fprintf('Nodes with thickness<ThickMin: ') ; fprintf('%i ',I)  ; fprintf('\n')
            fprintf('                    thickness: ') ; fprintf('%g ',F1.h(I))  ; fprintf('\n')
            if CtrlVar.ResetThicknessToMinThickness
                F1.h(F1.h<CtrlVar.ThickMin)=CtrlVar.ThickMin;
                %fprintf(CtrlVar.fidlog,' Found %-i thickness values less than %-g. Min thickness is %-g.',numel(indh0),CtrlVar.ThickMin,min(h));
                fprintf(CtrlVar.fidlog,' Setting h1(h1<%-g)=%-g \n ',CtrlVar.ThickMin,CtrlVar.ThickMin) ;
            end
        end
    end
    
    
    RunInfo.Forward.uvhActiveSetIterations(CtrlVar.CurrentRunStepNumber)=iActiveSetIteration-1 ;
    RunInfo.Forward.uvhActiveSetCyclical(CtrlVar.CurrentRunStepNumber)=isActiveSetCyclical;
    RunInfo.Forward.uvhActiveSetConstraints(CtrlVar.CurrentRunStepNumber)=numel(BCs1.hPosNode);
    
end
