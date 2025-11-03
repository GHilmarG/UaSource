






function [UserVar,RunInfo,h1,l1,BCs1]=MassContinuityEquationNewtonRaphson(UserVar,RunInfo,CtrlVar,MUA,F0,F1,l1,BCs1)

%% Solves: 
%
% $$\rho \, \frac{\partial h}{\partial t} + \nabla \cdot ( \rho \, \mathbf{v} h )  = \rho \,  a(h)$$
%
% subject to 
%
% $$h(x,y,t)  h_{\min} 0$$
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
nargoutchk(5,5)




MLC=BCs2MLC(CtrlVar,MUA,BCs1);
L=MLC.hL ; Lrhs=MLC.hRhs ;
if nargin==7 || isempty(l1.h) || (numel(l1.h)~=numel(Lrhs))
    l1.h=Lrhs*0 ;
end
dl=l1.h*0 ;
dh=F1.h*0;
BCsError=0;

% make sure initial point is feasible
F1.h(BCs1.hFixedNode)=BCs1.hFixedValue;
F0.h(BCs1.hFixedNode)=BCs1.hFixedValue;



iteration=0 ; rWork=inf ; rForce=inf; CtrlVar.NRitmin=0 ; gamma=1; rRatio=1; r=Inf;
RunInfo.Forward.hConverged=false;   BackTrackInfo.iarm=NaN;BackTrackInfo.Converged=0;

while true


    % exit criterion:
    %
    % First check if iteration>1, and back-step small and rRatio small. There could be a good reason for this. For example we
    % have already reached the NR minimum and no further reduction is possible due to rounding effect or because the cost
    % function is scaled in some ways that make further reduction irrelevant for the solution.


    if rRatio<0.01 && rWork>1e-20 % OK I'm hard wiring in here an option of additional iteration
        % The argument being that we might most likely still be in the second-order convergence]
        % and not limited by numerical rounding errors. (Consider taking this out once
        % all is fine.)

        ResidualsCriteria=false; % If there was a significant reduction in last iteration, do continue.

    elseif iteration>1 ...
            && (gamma < max(CtrlVar.hExitBackTrackingStepLength,CtrlVar.BacktrackingGammaMin) ...
            || rRatio>0.99)

        % This exit criterion applies if in last iteration either back-step or the reduction was very small. This indicates
        % difficulties with convergence, but this might simply be due to rounding errors making a further reduction impossible.

        ResidualsCriteria=(rWork<CtrlVar.hDesiredWorkAndForceTolerances(1)  && rForce<CtrlVar.hDesiredWorkAndForceTolerances(2))...
            && (rWork<CtrlVar.hDesiredWorkOrForceTolerances(1)  || rForce<CtrlVar.hDesiredWorkOrForceTolerances(2))...
            && iteration >= CtrlVar.NRitmin;


    else

        % This is the otherwise default exit criterion
        ResidualsCriteria=(rWork<CtrlVar.hAcceptableWorkAndForceTolerances(1)  && rForce<CtrlVar.hAcceptableWorkAndForceTolerances(2))...
            && (rWork<CtrlVar.hAcceptableWorkOrForceTolerances(1)  || rForce<CtrlVar.hAcceptableWorkOrForceTolerances(2))...
            && iteration >= CtrlVar.NRitmin;


    end


    if iteration>=CtrlVar.hSolverMaxIterations && (rRatio>0.9)
        if CtrlVar.InfoLevelNonLinIt>=1
            fprintf('hEquationNewtonRaphson: Maximum number of NR iterations (%i) reached. \n ',CtrlVar.hSolverMaxIterations)
        end
        RunInfo.Forward.hConverged=false;
        break
    end


    if ResidualsCriteria
        if CtrlVar.InfoLevelNonLinIt>=1
            fprintf('hEquationNewtonRaphson: NR iteration converged in %i iterations with rForce=%g and rWork=%g \n',iteration,rForce,rWork)
        end
        RunInfo.Forward.hConverged=true;
        break
    end


    if RunInfo.BackTrack.Converged==0 || gamma<CtrlVar.hExitBackTrackingStepLength
        if CtrlVar.InfoLevelNonLinIt>=1
            fprintf('hEquationNewtonRaphson: Backtracking stagnated or step too small. Exiting after %i iterations with rWork=%g and rForce=%g \n',iteration,rWork,rForce)
        end
        RunInfo.Forward.hConverged=false;
        break
    end

    iteration=iteration+1 ;



    [UserVar,R,K]=MassContinuityEquationAssembly(UserVar,RunInfo,CtrlVar,MUA,F0,F1) ;



    if ~isempty(L)

        frhs=-R-L'*l1.h        ;  % This needs to be identical to what is defined in the CalcCostFunctionhEquation
        grhs=Lrhs-L*F1.h; % Here the argument is that frhs has the units: [\varphi] area/time
        % while grhs has the units [\varphi], where [\varphi] are the units of
        % the level-set function itself, ie [\varphi]
    else
        frhs=-R;
        grhs=[];
    end

    [dh,dl]=solveKApe(K,L,frhs,grhs,dh,dl,CtrlVar);
    dh=full(dh);


    if any(isnan(dh))
        save TestSave
        error('hEquationNewtonRaphson:NaNinSolution','NaN in the solution for dh')
    end

    Func=@(gamma) CalcCostFunctionhEquation(UserVar,RunInfo,CtrlVar,MUA,gamma,F1,F0,L,Lrhs,l1.h,dh,dl);


    gamma=0 ; [r0,~,~,rForce0,rWork0,D20]=Func(gamma);
    gamma=1 ; r1=Func(gamma);

    slope0=-2*r0 ;

    if CtrlVar.hInfoLevel>=1000
        CtrlVar.InfoLevelBackTrack=1000 ;
    end


    CtrlVar.BacktrackingGammaMin=CtrlVar.hExitBackTrackingStepLength;
    nOut=6;
    varargout=cell(nOut-1,1);
    [gamma,r,BackTrackInfo,varargout{:}]=BackTracking(slope0,1,r0,r1,Func,CtrlVar,nOut);

    rForce=varargout{3};
    rWork=varargout{4};
    D2=varargout{5};
    %[r1Test,~,~,rForce,rWork,D2]=Func(gamma);



    if CtrlVar.hInfoLevel>=100 && CtrlVar.doplots==1
        nnn=30;
        gammaTestVector=zeros(nnn,1) ; rForceTestvector=zeros(nnn,1);  rWorkTestvector=zeros(nnn,1); rD2Testvector=zeros(nnn,1);
        Upper=2.2;
        Lower=-0.5 ;
        if gamma>0.7*Upper ; Upper=2*gamma; end
        parfevalOnAll(gcp(), @warning, 0, 'off','MATLAB:decomposition:genericError');
        for I=1:nnn
            gammaTest=(Upper-Lower)*(I-1)/(nnn-1)+Lower;
            [rTest,~,~,rForceTest,rWorkTest,D2Test]=Func(gammaTest);
            gammaTestVector(I)=gammaTest ; rForceTestvector(I)=rForceTest; rWorkTestvector(I)=rWorkTest;  rD2Testvector(I)=D2Test;
        end

        gammaZero=min(abs(gammaTestVector)) ;
        if gammaZero~=0
            [rTest,~,~,rForceTest,rWorkTest,D2Test]=Func(0);
            gammaTestVector(nnn+1)=0 ; rForceTestvector(nnn+1)=rForceTest; rWorkTestvector(nnn+1)=rWorkTest;  rD2Testvector(nnn+1)=D2Test;
        end

        [gammaTestVector,ind]=unique(gammaTestVector) ; rForceTestvector=rForceTestvector(ind) ; rWorkTestvector=rWorkTestvector(ind) ; rD2Testvector=rD2Testvector(ind) ;
        [gammaTestVector,ind]=sort(gammaTestVector) ; rForceTestvector=rForceTestvector(ind) ; rWorkTestvector=rWorkTestvector(ind) ; rD2Testvector=rD2Testvector(ind) ;


        SlopeForce=-2*rForce0;
        SlopeWork=-2*rWork0;
        SlopeD2=-D20;
        CtrlVar.MinimisationQuantity=CtrlVar.hMinimisationQuantity;
        [ForceFig,WorkFig]=PlotCostFunctionsVersusGamma(CtrlVar,RunInfo,gamma,r,iteration,"-h-",...
            gammaTestVector,rForceTestvector,rWorkTestvector,rD2Testvector,...
            SlopeForce,SlopeWork,SlopeD2,rForce,rWork,D2);

    end



    F1.h=F1.h+gamma*dh;
    l1.h=l1.h+gamma*dl;

    rRatio=r/r0;
    if CtrlVar.hInfoLevel>=1
        if ~isempty(L)
            BCsError=norm(Lrhs-L*F1.h);
        end

        if CtrlVar.InfoLevelNonLinIt>= 1
            fprintf(CtrlVar.fidlog,'NR-h:%3u/%-2u g=%-14.7g , r/r0=%-14.7g ,  r0=%-14.7g , r=%-14.7g , rForce=%-14.7g , rWork=%-14.7g , BCsError=%-14.7g \n ',...
                iteration,BackTrackInfo.iarm,gamma,rRatio,r0,r,rForce,rWork,BCsError);
        end
    end

    % Make sure to update s and b as well, but do not reset thickness within the non-linear loop

    CtrlVar.ResetThicknessToMinThickness=0;
    [F1.b,F1.s]=Calc_bs_From_hBS(CtrlVar,MUA,F1.h,F1.S,F1.B,F1.rho,F1.rhow);
    %[UserVar,F1]=GetMassBalance(UserVar,CtrlVar,MUA,F1); % actually this call only needed if mass-balance depends on h
end

h1=F1.h ; % Because I don't return F1


if numel(RunInfo.Forward.hIterations) < CtrlVar.CurrentRunStepNumber
    RunInfo.Forward.hIterations=[RunInfo.Forward.hIterations;RunInfo.Forward.hIterations+NaN];
    RunInfo.Forward.hResidual=[RunInfo.Forward.hResidual;RunInfo.Forward.hResidual+NaN];
    RunInfo.Forward.hBackTrackSteps=[RunInfo.Forward.hBackTrackSteps;RunInfo.Forward.hBackTrackSteps+NaN];
end

RunInfo.Forward.hIterations(CtrlVar.CurrentRunStepNumber)=iteration;
RunInfo.Forward.hResidual(CtrlVar.CurrentRunStepNumber)=r;
RunInfo.Forward.hBackTrackSteps(CtrlVar.CurrentRunStepNumber)=BackTrackInfo.iarm ;




end


