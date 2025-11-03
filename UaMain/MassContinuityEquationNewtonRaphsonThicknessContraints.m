




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
%% Solves: 
%
% $$\rho \, \frac{\partial h}{\partial t} + \nabla \cdot ( \rho \, \mathbf{v} h )  = \rho \,  a(h)$$
%
% subject to 
%
% $$h(x,y,t) > h_{\min}$$
%
% for $h$, using an implicit approach with respect to $h$.
%
% We can also write this as
%
% $$ \mathcal{M}(h) = 0 $$
%   
% where
%
% $$ \mathcal{M}(h):= \rho \, \frac{\partial h}{\partial t} + \nabla \cdot ( \rho \, \mathbf{v} h )  - \rho \,  a(h)$$
%
%
% The Galerkin weak form can then we written as
%
% $$ 0= \langle \mathcal{M}(h) \vert \phi_p \rangle = \int \mathcal{M}(h) \, \phi_p \; \mathrm{d} \mathcal{A} $$
%
% where all fields are expanded as
%
% for $$p=1 \ldots n $$, that is
%
%
% $$h(x,y) = \sum_{p=1}^n h_p \, \phi_p(x,y) $$
% 
% Here however the streamline-upwind Petrov-Galerkin method is used where 
%
% $$ 0= \langle \mathcal{M}(h) \vert \phi_p + \tau \, \mathbf{v} \cdot \nabla \phi_p  \rangle $$
%
% Therefore
%
% $$ 0= \left  \langle  \rho \, \partial_t h  + \nabla \cdot ( \rho \, \mathbf{v} h )  - \rho \,  a(h)  \; |  \; \phi_p + \tau \, \mathbf{v} \cdot \nabla \phi_p  \right \rangle $$
%
% As can be seen, the SUPG-perturbation is applied to all terms. This is therefore a consistent formulation. 
%
% Time discretization for $t=t_0$ to $t=t_1$ is done using the $\Theta$-method, i.e
% 
% $$\frac{\partial h}{\partial t} \approx \frac{1}{\Delta t} (h_1- h_0 ) $$
%
% with
%
% $$\Delta t= t_1-t_0 $$
%
% giving 
%
% $$ \rho \, \frac{1}{\Delta t} (h_1- h_0 ) = (1-\Theta) \left (  \rho\, a(h_0) - \nabla \cdot ( \rho \, \mathbf{v}_0 h_0 ) \right ) + \Theta \left (  \rho\, a(h_1) - \nabla \cdot ( \rho \, \mathbf{v}_1 h_1 ) \right )  $$
%
% After discretization the non-linear forward model can be written on the form
%
%
% $$ \mathbf{f}(\mathbf{h}) = \mathbf{0}  $$ 
%
% where
%
% $$ \mathbf{f} \in R^{n}$$
%
% and
%
% $$ \mathbf{h} \in R^n$$
%
% subject to the multi-nodal linear constraints
%
% $$
% \mathbf{L} \mathbf{h} = \mathbf{b}
% $$
%
% where  $$ \mathbf{L} \in R^{p\times n}$$ and  $$ \mathbf{b} \in R^{p}$$
%
% and furthermore subject to the positive thickness constraints
%
% $$ h_i > h_{\min} $$ for $$i=1 \ldots n$$
%
% 
%
% The multi-nodal (linear) constraints are enforced using the method of Lagrange multipliers. The positive thickness
% constraints can be enforced using the active set method and/or by applying a penalty term whereby an additional implicit
% mass balance term is added whenever $$a(h)< 0 $$. The resulting additional thickness constraints are internally added to
% the matrix $$\mathbf{L}$$, and updated after each Newton-Raphson solve.
%
% The Newton-Raphson system is then: 
%
% $$
% \left [ \begin{array}{cc}
% \mathbf{K} & \mathbf{L}^T  \\
% \mathbf{L} & \mathbf{0} 
% \end{array} \right ]
% \left [ \begin{array}{c}
% \Delta \mathbf{h} \\
% \Delta \mathbf{\lambda}
% \end{array} \right ]
% =\left [ \begin{array}{c}
% -\mathbf{f}(\mathbf{h}_k) - \mathbf{L}^T \mathbf{\lambda}_k \\
%   \mathbf{b} - \mathbf{L} \mathbf{h}_k
% \end{array} \right ]
% $$
%
% where $$k$$ is an iteration number (i.e. refers to the vectors at each iteration, and not to the individual components of
% the respective vectors).
%
% Here
%
% $$ \mathbf{K} = \mathrm{d} \mathbf{f}/\mathrm{d}\mathbf{h} $$
%
% and therefore
%
% $$ \mathbf{K} \in R^{n\times n}$$
%
% The Lagrange multipliers are
%
% $$ \mathbf{\lambda} \in R^{p}$$
%
% The Newton update is:
%
% $$ \mathbf{h}_{k+1} = \mathbf{h}_{k} + \Delta \mathbf{h} $$
%
% and
%
% $$ \mathbf{\lambda}_{k+1} = \mathbf{\lambda}_{k} + \Delta \mathbf{\lambda} $$
%
% where $k$ is the Newton-Raphson iteration step number.
%
% The iteration is continued until 
%
% $$
% \| \mathbf{f}(\mathbf{h}_k) + \mathbf{L}^T \mathbf{\lambda}_k \| < \epsilon
% $$
%
% Because the nodal constraints are linear, they are always fulfilled exactly at each iteration step (provided a full Newton
% step is taken).
%
%
%
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






