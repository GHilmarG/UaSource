



function [UserVar,RunInfo,F1,l1,BCs1,dt]=uvh(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1)


narginchk(9,9)


dt=CtrlVar.dt;


RunInfo.Forward.ActiveSetConverged=1;
RunInfo.Forward.uvhIterationsTotal=0;
iActiveSetIteration=0;
isActiveSetCyclical=NaN;
nlIt=nan(CtrlVar.NRitmax,1);


%%
% if any(F1.h<CtrlVar.ThickMin)
% 
%     indh0=find(F1.h<CtrlVar.ThickMin);
%     fprintf('uv: Found too small ice thicknesses in a prognostic forward run.\n')
%     fprintf('In total %-i negative ice thickness values found, with min ice thickness of %f. \n ',numel(indh0),min(F1.h));
% 
%     if CtrlVar.ResetThicknessToMinThickness==0
%         CtrlVar.ResetThicknessToMinThickness=1;
%     end
% 
%     fprintf('For the purpose of the uvh solve, these thickness values will be set to %f \n',CtrlVar.ThickMin)
%     [F1.b,F1.s,F1.h,F1.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F1.h,F1.S,F1.B,F1.rho,F1.rhow);
% 
% end
%%



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

        [UserVar,RunInfo,F1,F0,~,l1,BCs1,dt]=uvh2NotConvergent(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1) ;

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

    [UserVar,RunInfo,F1,l1,BCs1,~,Activated,Released]=ActiveSetInitialisation(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1) ;
    LastReleased=Released;  % Maybe here I should consider the possibility that the mesh has not changed and that I can use again the previous LastReleased and LastActivated list?
    LastActivated=Activated;
    iCounter=1;

    while true   % active-set loop

        fprintf(CtrlVar.fidlog,' ----------> Active-set Iteration #%-i <----------  \n',iActiveSetIteration);

        if CtrlVar.ThicknessConstraintsInfoLevel>=1
            if numel(BCs1.hPosNode)>0 || CtrlVar.ThicknessConstraintsInfoLevel>=10
                fprintf(CtrlVar.fidlog,' Number of active thickness constraints is %-i \n',numel(BCs1.hPosNode));
            end
        end

        % After some thought decided not to make iterate feasible, as the BCs should take care of this.
        % The issue is that by modifying h,ub,vb,  a repeat of a fully converged uvh solve does not start at
        % previous converged solution. Furthermore, with respect to the h thickness, the initial active set will not include those
        % nodes where the h has been modified, prior to solve, in this way. Most likely h will then again become negative at those
        % nodes and this will then likely result in an additional active set update.
        %
        %
        F1.h(BCs1.hPosNode)=CtrlVar.ThickMin;       % Make sure iterate is feasible, but this should be dealt with by the BCs anyhow
        F1.ub(BCs1.hPosNode)=F0.ub(BCs1.hPosNode);  % might help with convergence
        F1.vb(BCs1.hPosNode)=F0.vb(BCs1.hPosNode);

        [UserVar,RunInfo,F1,l1,BCs1]=uvh2D(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);

        if ~RunInfo.Forward.uvhConverged
            [UserVar,RunInfo,F1,F0,l0,l1,BCs1,dt]=uvh2NotConvergent(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l0,l1,BCs1);
        end

        if any(F1.h(BCs1.hPosNode) < CtrlVar.ThickMin)

            % Due to numerical errors it can sometimes happen that
            % on return even nodes in the active set
            % violate the pos thickness constraint.

            F1.h(BCs1.hPosNode) = CtrlVar.ThickMin ;
            if CtrlVar.ThicknessConstraintsInfoLevel>=1
                fprintf('Active-set: Resetting to limit. \n')
            end
        end




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

        % Update active set once uvh solution has been found.
        % Note: If this active-set iteration does not converge,
        %       the uvh solution and the active set are not consistence,
        %       i.e. the uvh solution was not obtained for the BCs in the active set.
        [UserVar,RunInfo,BCs1,~,isActiveSetModified,isActiveSetCyclical,Activated,Released]=ActiveSetUpdate(UserVar,RunInfo,CtrlVar,MUA,F1,l1,BCs1,iActiveSetIteration,LastReleased,LastActivated);

        LastReleased=Released;
        LastActivated=Activated;
        iActiveSetIteration=iActiveSetIteration+1;

        if ~isActiveSetModified
            fprintf(' Leaving active-set loop because active set unchanged in last active-set iteration. \n')
            break

        end

        if isActiveSetCyclical

            fprintf(" Leaving active-set pos. thickness loop because it has become cyclical\n");

            % ToDo (1 Sept 2024):  Do I want to do one final uvh solve since the active set has changed?

            % [UserVar,RunInfo,F1,l1,BCs1]=uvh2D(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);


            break

        end


        if  iActiveSetIteration > CtrlVar.ThicknessConstraintsItMax
            RunInfo.Forward.ActiveSetConverged=0;
            fprintf(' Leaving active-set pos. thickness loop because number of active-set iteration (%i) greater than maximum allowed (CtrlVar.ThicknessConstraintsItMax=%i). \n ',iActiveSetIteration,CtrlVar.ThicknessConstraintsItMax)
            break
        end

    end   % end of active set loop


    % nodes where thickness still is too small
    ToBeActivated=setdiff(find(F1.h < CtrlVar.ThickMin),BCs1.hFixedNode) ;  % This should not really be needed as this must be equal to the "Activated" set, except that Active set is not updated if
    % number of to-be-activated nodes less than CtrlVar.MinNumberOfNewlyIntroducedActiveThicknessConstraints


    if numel(ToBeActivated)>0

        warning('some h1 <ThickMin on return from active-set loop. min(h1)=%-g',min(F1.h)) ;

        if CtrlVar.ThicknessConstraintsInfoLevel >= 10
            fprintf('\tNodes with thickness<ThickMin: ') ; fprintf('\t %10i \t %10i \t %10i \t %10i \t  %10i \t  %10i \t  %10i \t  %10i \t  %10i \t  %10i \n  ',ToBeActivated)  ; fprintf('\n')
            fprintf('\t                    thickness: ') ; fprintf('\t %10g \t %10g \t %10g \t %10g \t  %10g \t  %10g \t  %10g \t  %10g \t  %10g \t  %10g \n  ',F1.h(ToBeActivated))  ;  fprintf('\n')
        end

        if CtrlVar.ResetThicknessToMinThickness
            F1.h(F1.h<CtrlVar.ThickMin)=CtrlVar.ThickMin;
            %fprintf(CtrlVar.fidlog,' Found %-i thickness values less than %-g. Min thickness is %-g.',numel(indh0),CtrlVar.ThickMin,min(h));
            fprintf(CtrlVar.fidlog,' Setting h1(h1<%-g)=%-g \n ',CtrlVar.ThickMin,CtrlVar.ThickMin) ;
        end

    end


    if  ~RunInfo.Forward.ActiveSetConverged
        fprintf(CtrlVar.fidlog,' Warning: In enforcing thickness constraints and finding a critical point, the loop was exited due to maximum number of iterations (%i) being reached. \n',CtrlVar.ThicknessConstraintsItMax);
    else
        if numel(BCs1.hPosNode)>0
            fprintf(CtrlVar.fidlog,'----- Active-set iteration converged after %-i iterations, constraining %-i thicknesses  \n \n ',iActiveSetIteration-1,numel(BCs1.hPosNode));
        end
    end



end



RunInfo.Forward.uvhActiveSetIterations(CtrlVar.CurrentRunStepNumber)=iActiveSetIteration-1 ;
RunInfo.Forward.uvhActiveSetCyclical(CtrlVar.CurrentRunStepNumber)=isActiveSetCyclical;
RunInfo.Forward.uvhActiveSetConstraints(CtrlVar.CurrentRunStepNumber)=numel(BCs1.hPosNode);

F1.solution="-uvh-" ; 

end
