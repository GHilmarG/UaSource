





function [I,dIdp,Psi_x,Psi_y,F,dFduv]=Misfit(CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,dFduv)

narginchk(9,9)
nargoutchk(6,6)


%%
%
%           J(q(p),p) = I(q(p)) + R(p)
%
% Forward model:   f(q(p),p)=0
%
%
% Calculate data misfit functions and gradients with respect to control
% variable p containing and the state variable q.
%
% The control variable p can involve A, B, C.
%
% The adjoint method is used to calculate the derivatives
%
% The misfit term, $I$, on the form
%
% $$I= \|u-u_{\mathrm{Meas}} \| +  \|v-v_{\mathrm{Meas}} \| + \| \dot{h}-\dot{h}_{\mathrm{Meas}} \| $$
%
% where
%
% $$ \dot{h} = a - \partial_x (u h ) - \partial_y (v h) $$
%
% that is
%
% $$I= \|u-u_{\mathrm{Meas}} \| +  \|v-v_{\mathrm{Meas}} \| + \| (\left ( a - \partial_x (u h ) - \partial_y (v h)  \right ) -\dot{h}_{\mathrm{Meas}} \| $$
%
% Thus, the misfit term is not an explicit function of $\dot{h}$
%
% $$I=I(u(p),v(p)) $$
%
% where $p$ are the parameters we want to invert for, i.e. any or all of $A$, $B$ and $C$.
%
% Since $\dot{h}$ is already expressed directly in terms of $u$ and $v$ using the mass-conservation equation, the forward
% model consists of the momentum equations only.
%
% Solving the forward model is
%
% $$\partial F /\partial u = \partial I/\partial u $$
%
% The partial derivatives of $I$ are straightforward to calculate, apart from
%
% $$\partial I_{\dot{h}} / \partial u$$
%
% which more correctly should be written as
%
% $$\delta I_{\dot{h}} \, \delta u $$
%
%%




DAI=[];
DBI=[];
DCI=[];
dIdp=[] ;

Psi_x=[]; Psi_y=[]; 


[is_uv_meas,is_dhdt_meas]=is_uv_dhdt_Meas(CtrlVar);
[isA,isB,isC] = isABC(CtrlVar) ;

if is_dhdt_meas
    [~,F.dhdt]=dhdtExplicit([],CtrlVar,MUA,F,BCs);
end


uErr=full(sqrt(spdiags(Meas.usCov)));
usres=(F.ub-Meas.us)./uErr;
vErr=full(sqrt(spdiags(Meas.vsCov)));
vsres=(F.vb-Meas.vs)./vErr;

%% Calculate misfit term I and its (explicit) gradients with respect to the state variable q, i.e. u and v
%
% This is straightforward as the data misfit terms are explicit functions of the state variables u and v.
%


Iuv=0 ; Ihdot=0 ;
duIdu=zeros(MUA.Nnodes,1) ;
dvIdv=zeros(MUA.Nnodes,1) ;

duIhdot=zeros(MUA.Nnodes,1);
dvIhdot=zeros(MUA.Nnodes,1);
dhIhdot=zeros(MUA.Nnodes,1) ;


if ~isfield(MUA,"M")
    MUA.M=MassMatrix2D1dof(MUA);
end

if is_uv_meas

    duIdu=(MUA.M*usres)./uErr/MUA.Area;    %    usres=(us-Meas.us)./uErr;
    dvIdv=(MUA.M*vsres)./vErr/MUA.Area;
    Iuv=full(usres'*MUA.M*usres+vsres'*MUA.M*vsres)/2/MUA.Area;

end

% Derivatives of the Ihdot misfit term. This is a bit more complicated as I need to take the derivative of the
% mass-conservation equation with respect to v.
%
% BTW, as the regularization term (R), does not depend on v these I derivatives with respect to v are the derivatives of J
% with respect to v.
%
% There is an interesting added aspect to this, because the mass-conservation involves both v and B, i.e. both the output of
% the forward model (v), and the variable which we are inverting for (B), the 'output' of the inverse model. For this reason,
% we have here a contribution from the misfit term (I) to the quantity <\partial_B J, \phi > which otherwise would only arise
% from the regularization term (R).  Therefore, we must here calculate the derivative of the cost function term with respect
% to both u and B. This addition is only needed when inverting for B while also including the hdot cost function term.
%
if is_dhdt_meas
    % when including dhdt meas, there is an extra contribution to the RHS of the adjoint system related to
    %
    % $$ \delta_u J_{\dot{h}} $$ and  $$ \delta_v J_{\dot{h}} $$ and
    %
    %
    [Ihdot,duIhdot,dvIhdot,dhIhdot]=EvaluateJhdotAndDerivatives([],CtrlVar,MUA,F,BCs,Meas);

end

I=Iuv+Ihdot ;  %
duvIduv=[duIdu(:)+duIhdot(:);dvIdv(:)+dvIhdot(:)];

if CtrlVar.TestAdjointFiniteDifferenceType=="complex step differentiation"
    CtrlVar.TestForRealValues=false;
end



%% Calculate the (implicit) gradients with respect to the control variable q, i.e. A, B and C
%
%
% Calculate the gradient of the misfit function I with respect to the control variables (model parameters) p (here A and B or C).
%
% This is a bit tricky because I=I(u(p))
%
if CtrlVar.Inverse.CalcGradI


    switch lower(CtrlVar.Inverse.DataMisfit.GradientCalculation)

        case {"fixpoint","fixpointc","-fixpoint-","-fixpointc-"}

            error("Misfit:OptionNoLongerSupported","fixpoint gradients now obsolete and no longer supported")


        case {"adjoint","-adjoint-"}


            %% Inverse problem
            %
            % Forward model:
            %   f(u(p),p)=0
            %
            %% Step 1: solve the non-linear forward problem
            %

            % We need a fully converged solution here
            % This will have been done ahead of the call, but I found that doing this (possibly again) increases accuracy.
            % Anyhow, it might be the case that just the Misfit.m is called from somewhere else, so one should really do this here as
            % well.
            [~,~,F,l,dFduv]= uv([],[],CtrlVar,MUA,BCs,F,l);

            %% Step 2:  Solve adjoint equation, i.e.   dfuv l = -dJduv
            %
            % $$ \langle \Psi \vert \delta_q \mathcal{F}[\phi_k] \rangle =- \delta_q J[\phi_k] $$
            %
            % fprintf(' Solve adjoint problem \n ')
            % I need to impose boundary conditions on lx and ly
            % if the problem is (fully) adjoint I have exactly the same BC
            % I need to solve
            %
            % [Kxu Kxv Luv'] [lx]        =  [ u-uMeas ]
            % [Kyu Kyv     ] [ly]        =  [ v-vMeas ]
            % [  Luv      0] [lambdauv]     [ Luvrhs  ]
            % All matrices are Nnodes x Nnodes, apart from:
            % Luv is #uv constraints x 2 Nnodes

            MLC_Adjoint=BCs2MLC(CtrlVar,MUA,BCsAdjoint);
            LAdjoint=MLC_Adjoint.ubvbL;
            LAdjointrhs=MLC_Adjoint.ubvbRhs;
            lAdjoint=zeros(numel(LAdjointrhs),1) ;

            %duvJ=duvIduv;     % Because this is the only J term that depends on (u,v).
            RHS_Adjoint=-duvIduv;
            % If the regularization term also depended on the measurements q, ie R=R(u,v) then this would not be correct.

            % Now solve the linear adjoint problem for lambda
            [lambda,lAdjoint]=solveKApeSymmetric(dFduv,LAdjoint,RHS_Adjoint,LAdjointrhs,[],lAdjoint,CtrlVar);


            Psi_x=real(lambda(1:MUA.Nnodes)) ;
            Psi_y=real(lambda(MUA.Nnodes+1:2*MUA.Nnodes));

         


            %% Step 3:  <d_p F^* \lambda>,
            %
            % Note that I'm adding the d_p R term in the regularization step.
            %
            % But I need to include a possible <d_p I , \phi> term here.
            %
            % For p=A and p=C, d_p I =0 because I is not an explicit function of A and C
            %
            % But for b, d_b I = p_x (u db)

            if isC

                % $$ \langle  \delta_{C_i} F^x \phi_i | \Psi_x \rangle + \langle  \delta_{C_i} F^y \phi_i| \Psi_y \rangle $$
                dCFuvLambda=dIdCq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);

                dCI=0 ;      % This is the explicit derivative of the misfit term, I, with respect to C. There is no such dependency here
                % as the misfit term I is not an explicit function of C, so this equals to zero.
                % Note, that there is an explicit dependency on C in the regularization term, but this is added elsewhere (in the
                % Regularisation.m function)

                DCI=dCFuvLambda+dCI;  % this is the part of the dI/dC derivative which is due to the implicit dependency
                % of I on C because the velocities depend on C,

            end

            if isA


                dAFuvLambda=dIdAq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y);

                dAI=0 ; % No explicit dependency of the misfit term I on A.

                DAI=dAFuvLambda+dAI;
            end


            if isB

                OnlyGrounded=true;


                if OnlyGrounded
                    %  p= B ;

                    dBdp=  1+zeros(MUA.Nnodes,1);
                    dbdp=  F.GF.node ; % - (1-F.GF.node).*F.GF.node.*F.rho/F.rhow;
                    dhdp= -F.GF.node ;

                else

                    %  p= B ;
                    dBdp=  F.GF.node ; %
                    dbdp=  F.GF.node ; % - (1-F.GF.node).*F.GF.node.*F.rho/F.rhow;
                    dhdp= -F.GF.node ;



                end

                % dIdB= dhF^* \lambda + dhJ
                % if only -dhdt- meas and no regularization
                % then dJdB=dh/db*dhJhdot


                dBFuvLambda=dIdbq(CtrlVar,MUA,F,BCs,BCsAdjoint,Psi_x,Psi_y,dhdp,dbdp,dBdp);
                %   dBFuvLambda2=dIdBq2(CtrlVar,MUA,uAdjoint,vAdjoint,F);
                %dBFuvLambda=dBFuvLambda2;


                dBI=dhdp.*dhIhdot;  % The Ihdot misfit term includes an explicit dependency on B, which is here accounted for.

                dBI=ApplyAdjointGradientPreMultiplier(CtrlVar,MUA,BCsAdjoint,CtrlVar.Inverse.AdjointGradient.UseBCs.B,dBI); % added 7 Jan 2025

                DBI=dBFuvLambda+dBI;
                %
                % if CtrlVar.Inverse.OnlyModifyBedUpstreamOfGL
                %     F.GF=IceSheetIceShelves(CtrlVar,MUA,F.GF,GLgeo,GLnodes,GLele) ;
                %     DBI(~F.GF.NodesUpstreamOfGroundingLines)=0;
                % end
                %

                % UaPlots(CtrlVar,MUA,F,dBFuvLambda,FigureTitle="dBFuvLambda")
                % UaPlots(CtrlVar,MUA,F,dBI,FigureTitle="dBI") ;  CM=cmocean('balanced',25,'pivot',0) ; colormap(CM);

            end



        otherwise

            fprintf(" CtrlVar.Inverse.DataMisfit.GradientCalculation has the value %s \n",CtrlVar.Inverse.DataMisfit.GradientCalculation)
            fprintf(" but the only allowed values are ''fixpoint'' and ''adjoint'' \n")
            error(" which case? ")

    end

    dIdp=[DAI;DBI;DCI] ;  % 2026 Feb
end




end


