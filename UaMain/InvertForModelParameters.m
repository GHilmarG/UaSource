




function [UserVar,F,l,InvFinalValues,RunInfo]=...
    InvertForModelParameters(UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo)

narginchk(11,11)

%%
%
%
% As shown in detail in the Ua compendium, an inversion for the parameters p, can be formulated as a constrained minimization
% problem:
% 
% $$J(p)=I(q(p),p) + R(p) $$
% 
% subject to
%
% $$F(q(p),p)=0$$
% 
% To use a gradient-based optimization method, we need to be able to calculate the derivatives of J with respect to p.
%
% Using the Lagrange method, we form the extended functional
%
% $$\mathcal{L}(p)=I(q(p),p) + R(p) + \langle F(q(p),p) | \lambda \rangle $$
%
%
% We can then find the directional derivative of $J(p)$ with respect to $p$ in the direction $\phi$ as follows:
%
% 1) For some $p$, solve for $q$ using the forward model:
%
% $$F(q(p),p)=0$$
%
% 2) Then for this $q$, solve the linear adjoint problem:
%
% $$ \langle (d_q F)^* \lambda | \phi \rangle = - \langle  d_q J | \phi \rangle  $$
%
% 3) And then evaluate the total derivative with respect to $p$ as
%
% $$ d_p J = \langle (d_q F)^* | \lambda  \rangle + \partial_p J $$
%
% As an example, consider a $B$ inversion using momentum and mass conservation (for grounded ice where $h=s-B$):
%
% 
% $$J(B)= I(v(B),B) + R(B) $$
%
% subject to
%
% $$P(v(B),B) = 0 $$   (momentum)
%
% $$M(B) = 0 $$   (mass)
% 
% we could introduce Lagrange multipliers for both of these equations, but we can also use the fact that $\dot{h}$ can so
% easily be calculated from the mass conservation equation, and put the mass-conservation directly into the cost function $J$
%
% $$J(B)= \| v_c - v_m \| + \| \dot{h}_c - \dot{h}_m \|  + \|B_c - B_m \| + \langle P(v(B),B) | \lambda \rangle $$
%
% we can write this as
%
% $$J(B)= I_v + I_{\dot{h}}  + R(B) + \langle P(v(B),B) | \lambda \rangle $$
%
% were $I_v$ and $I_{\dot{h}}$ are misfit terms,  and $R$ a regularization term, and where we simply calculate evaluate
% 
% $$\dot{h}_c=a - \nabla (v \, (s-B) ) $$
%
% and insert into the misfit term, i.e.
%
% $$I_{\dot{h}} = \| (a - \nabla (v \, (s-B) ) - \dot{h}_m \| $$
%
% The misfit term $I_{\dot{h}}$ is an explicit function of $v$ (i.e. the q variable).
%
% The (linear) adjoint problem now reads
%
% $$ \langle (d_v F)^* \lambda | \phi \rangle = - \langle  d_v I_v | \phi \rangle  - \langle  d_v I_{\dot{h}} | \phi \rangle
% $$
%
% The right-hand term is the derivative of $J=I+R$ with respect to $v$, but as $R$ does not depend on $v$, we only have those
% two terms involving $I$.
%
% The directional derivative of $J$ with respect to $B$ is calculated from
%
%
% $$ d_B J = \langle (d_B F)^* | \lambda \rangle + \partial_B J $$
%
%%

if isempty(CtrlVar.Inverse.InitialLineSearchStepSize) ||  CtrlVar.Inverse.InitialLineSearchStepSize==0
    CtrlVar.Inverse.InitialLineSearchStepSize=InvStartValues.SearchStepSize;
end

if CtrlVar.Inverse.Regularize.Field=="-logAGlen-logC-"

    if isempty(CtrlVar.Inverse.Regularize.logAGlen.ga)

        fprintf("The variable CtrlVar.Inverse.Regularize.logAGlen.ga is undefined! This variable needs to be defined in DefineInitialInputs.m \n")
        error("Input variable not defined.")
    end

    if isempty(CtrlVar.Inverse.Regularize.logAGlen.gs)

        fprintf("The variable CtrlVar.Inverse.Regularize.logAGlen.gs is undefined! This variable needs to be defined in DefineInitialInputs.m \n")
        error("Input variable not defined.")
    end

end



%% Define inverse parameters and anonymous function returning objective function, directional derivative, and Hessian
%

% F should always be populated with the fields used at each stage of the inversion, and all the calculations are done using
% F. Start by populating F with the starting values 
F=InvStartValues2F(CtrlVar,MUA,F,InvStartValues,Priors,Meas) ;
[F.b,F.s,F.h,F.GF]=Calc_bs_From_hBS(CtrlVar,MUA,F.h,F.S,F.B,F.rho,F.rhow);


F.GF=IceSheetIceShelves(CtrlVar,MUA,F.GF) ;

% p is the vector of the control variables, currently p=[A,B,C]
% with A, B or C here only being nonempty when inverted for, 
% This mapping between A, B and C into the control variable is done by F2p

% Make sure initial point is feasible
F.AGlen=kk_proj(F.AGlen,F.AGlenmax,F.AGlenmin) ;
F.C=kk_proj(F.C,F.Cmax,F.Cmin) ;
F.B=kk_proj(F.B,F.Bmax,F.Bmin) ;

% The parameter that we are inverting for are contained in the variable p. p0 is the starting value.
[p0,plb,pub]=F2p(CtrlVar,MUA,F); 


CtrlVar.Inverse.ResetPersistentVariables=1;


% JGH: Returns the cost function (J), the gradient of the cost function with respect to p (dJdp), and the Hessian (ddJddp).
% The Hessian of the regularization term (R) can usually be calculated exactly, while the Hessian of the misfit/likelihood term
% (I), can not. However, one can come up with an educated guess for the Hessian of I with respect to C.
[J0,dJdp,Hessian,JGHouts,F,RunInfo]=JGH(p0,plb,pub,UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);
CtrlVar.Inverse.ResetPersistentVariables=0;
% The parameters passed in the anonymous function are those that exist at the time the anonymous function is created.



% Function handles are created to the functions calculating the cost function, J, the gradient, dJdp, and the Hessian.
% This is then passed to the optimization libraries. For some reason the MATLAB optimization library requires a separate
% handle to the Hessian.

func=@(p) JGH(p,plb,pub,UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);   % returns the cost (J), gradient (G) and Hessian (H)
Hfunc=@(p,lambda) HessianAC(p,lambda,plb,pub,UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo); % returns the Hessian (H). 
                                                                                                                         % Somewhat annoyingly MATLAB optimisation toolbox 
                                                                                                                         % wants the Hessian returned in 
                                                                                                                         % a separate function, so I can't use JGH (!?).
                                                                                                                         % The function HessianAC is just a wrapper around 
                                                                                                                         % JGH and returns the same Hessian as JGH.

fprintf('\n +++++++++++ At start of inversion:  \t J=%-g \t I=%-g \t R=%-g  |grad|=%g \n \n',J0,JGHouts.MisfitOuts.I,JGHouts.RegOuts.R,norm(dJdp))

dJdpTest=[];

%%

if CtrlVar.Inverse.TestAdjoint.isTrue
    %% The correctness of the gradient calculation can be tested by comparing it with a brute-force finite differences calculations. 
   
    % Get the gradient using the adjoint method
    [J,dJdp,Hessian,JGHouts]=func(p0);
    
    
    %NA=numel(InvStartValues.AGlen);  % Number of parameters to invert for
    NA=MUA.Nnodes; 
    
    % Find the subset (iRange) in p, for which the brute-force gradient is to be calculated
    if isempty(CtrlVar.Inverse.TestAdjoint.iRange)
        iRange=1:NA;  % If iRange is left empty, do for all of p, i.e. with respect to values over all nodes 
                      
    else
        iRange=CtrlVar.Inverse.TestAdjoint.iRange;
    end
    
    % if the inversion is done for more than one field, then expand iRange accordingly. 
    switch strlength(CtrlVar.Inverse.InvertForField)
        
        case 2
            iRange=[iRange(:);iRange(:)+NA];
        case 3
            iRange=[iRange(:);iRange(:)+NA;iRange(:)+2*NA];
    end
    
    I=(iRange>=1) & (iRange <= numel(p0));  % As far as I can see, this should not be needed...
    iRange=iRange(I);
    
    % calculate brute force gradient

    % Gradient calculated using a brute-force finite difference approach 
    dJdpTest = CalcBruteForceGradient(func,p0,CtrlVar,iRange);

    filename=CtrlVar.Experiment+"BruteForceGradient";
    fprintf('BruteForceGradient save in the file : %s \n',filename)
    save(filename,'CtrlVar','UserVar','MUA','F','dJdpTest','iRange')
    
    
else
    
    
    %%
    
    if contains(CtrlVar.Inverse.MinimisationMethod,"Ua")
        
        [p,UserVar,RunInfo]=UaOptimisation(UserVar,CtrlVar,RunInfo,MUA,func,p0,plb,pub);
        
        
    elseif contains(CtrlVar.Inverse.MinimisationMethod,"Matlab")
        
        clear fminconOutputFunction fminconHessianFcn fminuncOutfun
        
        [p,RunInfo]=InversionUsingMatlabOptimizationToolbox3(UserVar,CtrlVar,RunInfo,MUA,func,p0,plb,pub,Hfunc);
        
        
    else
        
        fprintf(' CtrlVar.Inverse.MinimisationMethod has the value %s \n',CtrlVar.Inverse.MinimisationMethod)
        fprintf(' but can only have the values ''MatlabOptimization'' or ''UaOptimization''\n')
        error('what case? ')
    end
    
    % Here the final values from inversion, which are in the vector p, are copied across to the corresponding fields of F
    F=p2F(CtrlVar,MUA,p,F,Meas,Priors);
   
    % And a final additional call is made to get the cost function, J, and the gradient, dJdp, at the end of the optimization. In
    % principle, I guess it should be possible to get this information from the (external) optimization subroutine, but I don't
    % know how...
    [J,dJdp,Hessian,JGHouts,F,RunInfo]=JGH(p,plb,pub,UserVar,CtrlVar,MUA,BCs,F,l,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo);
    fprintf('\n +++++++++++ At end of inversion:  \t J=%-g \t I=%-g \t R=%-g  |grad|=%g \n \n',J,JGHouts.MisfitOuts.I,JGHouts.RegOuts.R,norm(dJdp))
    
    
end

% Put RAa, RAs, RCa, RCs in InvFinalValues
InvFinalValues=Vars2InvValues(CtrlVar,F,InvStartValues,J,dJdp,JGHouts,RunInfo,dJdpTest); 


end
