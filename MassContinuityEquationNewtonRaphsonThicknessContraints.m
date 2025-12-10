




function [UserVar,RunInfo,h1,l1,BCs1]=MassContinuityEquationNewtonRaphsonThicknessContraints(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1)

narginchk(8,8)
nargoutchk(5,5)


if CtrlVar.ThicknessConstraints



    iActiveSetIteration=0;

    % CtrlVar.ThicknessConstraintsItMax=5;
    % CtrlVar.ThicknessConstraintsInfoLevel=1;

    [UserVar,RunInfo,F1,l1,BCs1,isActiveSetModified,Activated,Released]=ActiveSetInitialisation(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1) ;
    LastReleased=Released;  % Maybe here I should consider the possibility that the mesh has not changed and that I can use again the previous LastReleased and LastActivated list?
    LastActivated=Activated;


    while true

        iActiveSetIteration=iActiveSetIteration+1;
        fprintf("h-Solve Nr. %i \n ",iActiveSetIteration)
        [UserVar,RunInfo,h1,l1,BCs1]=MassContinuityEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);



        F1.h=h1;


        fprintf("Active set update Nr. %i \n ",iActiveSetIteration)

        [UserVar,RunInfo,BCs1,lambdahpos,isActiveSetModified,isActiveSetCyclical,Activated,Released]=ActiveSetUpdate(UserVar,RunInfo,CtrlVar,MUA,F1,l1,BCs1,iActiveSetIteration,LastReleased,LastActivated);
        LastReleased=Released;
        LastActivated=Activated;




        if ~isActiveSetModified
            fprintf(' Leaving active-set loop because active set unchanged in last active-set iteration. \n')
            break

        end

        if isActiveSetCyclical

            fprintf(" Leaving active-set pos. thickness loop because it has become cyclical\n");

            fprintf("h-Solve Nr. %i \n ",iActiveSetIteration+1)
            [UserVar,RunInfo,h1,l1,BCs1]=MassContinuityEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1);

            break

        end


        if  iActiveSetIteration > CtrlVar.ThicknessConstraintsItMax
            RunInfo.Forward.ActiveSetConverged=0;
            fprintf(' Leaving active-set pos. thickness loop because number of active-set iteration (%i) greater than maximum allowed (CtrlVar.ThicknessConstraintsItMax=%i). \n ',iActiveSetIteration,CtrlVar.ThicknessConstraintsItMax)
            break
        end


    end




else


    [UserVar,RunInfo,h1,l1,BCs1]=MassContinuityEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1) ;


end



end






