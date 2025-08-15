




function [UserVar,RunInfo,h1,l1,BCs1]=MassContinuityEquationNewtonRaphsonThicknessContraints(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1)


%% Solves: 
%
% $$\rho \frac{\partial h}{\partial t} + \nabla \cdot ( \rho \, \mathbf{v} h )  =  \rho \, a(h)$$
%
% for $h$, using an implicit approach with respect to $h$.
%
% $\mathbf{v}$=(F.ub,F.vb)
%
% The approach use the active set to enforce $h>h_{\min}$, provided 
%
%   CtrlVar.ThicknessConstraints=true
% 
% If
%
%   CtrlVar.LevelSetMethod=true ; 
%
%   CtrlVar.LevelSetMethodAutomaticallyApplyMassBalanceFeedback=true;
%
% an additional mass-balance term is added, using same approach as done in the uvh solve.
%
% Also, an additional mass balance term can be added, as in the -uvh- solve, as a penalty term for too small ice thicknesses.
% This option is activated by setting:
%
%   CtrlVar.ThicknessPenalty=true ;
%
% Further details are provided in Ua2D_DefaultParameters.m
%
% The assembly is identical to the Khh assembly in the -uvh- solve, and does include gradients in density.
% 
%%


narginchk(8,8)
nargoutchk(3,5)

MUA=UpdateMUA(CtrlVar,MUA);  

if CtrlVar.ThicknessConstraints



    iActiveSetIteration=0;

  

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






