function [I,dIdp,ddIdpp,MisfitOuts]=Misfit(UserVar,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,RunInfo,dfuv)

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
%
%
%%


Area=MUA.Area;

DAI=[];
DBI=[];
DCI=[];

dIdp=[] ;
ddIdpp=sparse(1,1);

ddIdAA=[];
ddIdCC=[];


MisfitOuts.dIduv=[];
MisfitOuts.uAdjoint=[];
MisfitOuts.vAdjoint=[];

us=F.ub+F.ud ;
vs=F.vb+F.vd ;

%% Do some test on inputs, check if error covariance matrices are correctly defined
% calculate residual terms, i.e. difference between measured and modeled values.
if contains(CtrlVar.Inverse.Measurements,"-uv-","IgnoreCase",true)
    
    if isempty(Meas.us)
        fprintf('Meas.us is empty! \n')
        fprintf('Meas.us cannot be empty when inverting using surface velocities as data.\n')
        fprintf('Define Meas.us in DefineInputsForInverseRun.m \n')
        error('Misfit:us','Meas.us is empty.')
    end
    
    
    if isempty(Meas.vs)
        fprintf('Meas.vs is empty! \n')
        fprintf('Meas.vs cannot be empty when inverting using surface velocities as data.\n')
        fprintf('Define Meas.vs in DefineInputsForInverseRun.m \n')
        error('Misfit:us','Meas.vs is empty.')
    end
    
    if isdiag(Meas.usCov)
        uErr=sqrt(spdiags(Meas.usCov));
        usres=(us-Meas.us)./uErr;
    else
        error('Misfit:Cov','Data covariance matrices but be diagonal')
    end
    
    if isdiag(Meas.vsCov)
        vErr=sqrt(spdiags(Meas.vsCov));
        vsres=(vs-Meas.vs)./vErr;
    else
        error('Misfit:Cov','Data covariance matrices must be diagonal')
    end
end


if contains(CtrlVar.Inverse.Measurements,'-dhdt-','IgnoreCase',true)
    
    if isempty(Meas.dhdt)
        fprintf('Meas.dhdt is empty! \n')
        fprintf('Meas.dhdt cannot be empty when inverting using dhdt as data.\n')
        fprintf('Define Meas.dhdt in DefineInputsForInverseRun.m \n')
        error('Misfit:dhdt','Meas.dhdt is empty.')
    end
    
    
    if isempty(Meas.dhdtCov)
        fprintf('Meas.dhdtCov is empty! \n')
        fprintf('Meas.dhdtCov cannot be empty when inverting using dhdt as data.\n')
        fprintf('Define Meas.dhdtCov in DefineInputsForInverseRun.m \n')
        error('Misfit:dhdt','Meas.dhdt is empty.')
    end
    
    [UserVar,F.dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F,BCs);
    
    if ~isdiag(Meas.dhdtCov)
        error('Misfit:Cov','Data covariance matrices must be diagonal')
    end
    
end

%% Calculate misfit term I and its gradient with respect to the state variable q, i.e. u and v
%
% This is straightforward as the misfit term is an explicit function of u and v.
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

if contains(CtrlVar.Inverse.Measurements,"-uv-")
    
    duIdu=(MUA.M*usres)./uErr/Area;
    dvIdv=(MUA.M*vsres)./vErr/Area;
    Iuv=full(usres'*MUA.M*usres+vsres'*MUA.M*vsres)/2/Area;
    
end

% Derivatives of the Ihdot misfit term. This is a bit more complicated as I need to take the derivative of the
% mass-conservation equation with respect to v.
%
% BTW, as the regularization term (R), does not depend on v these I derivatives with respect to v are the derivatives of J
% with respect to v.
%
% There is an interesting added aspect to this, because the mass-conservation involves both v and B, i.e. both the output of
% the forward model (v), and the variable which we are inverting for (B), the 'output' of the inverse model. For this reason,
% we have here a contribution from the misfit term (I) to the quantity <\partial_B J, \phi > which otherwise would only arise from
% the regularization term (R).  Therefore, we must here calculate the derivative of the cost function term with respect to
% both u and B. This addition is only needed when inverting for B while also including the hdot cost function term.
%
if contains(CtrlVar.Inverse.Measurements,"-dhdt-")
    
    [Ihdot,duIhdot,dvIhdot,dhIhdot]=EvaluateJhdotAndDerivatives(UserVar,CtrlVar,MUA,F,BCs,Meas);
    
end


I=Iuv+Ihdot ;  %  
duvIduv=[duIdu(:)+duIhdot(:);dvIdv(:)+dvIhdot(:)];

if CtrlVar.TestAdjointFiniteDifferenceType=="complex step differentiation"
    CtrlVar.TestForRealValues=false;
end

if ~isreal(duvIduv)  && CtrlVar.TestForRealValues
    save TestSave ; error("MisfitFunction:dIduvNoReal","dIduv is not real! Possibly a problem with covariance of data.")
end


MisfitOuts.dIduv=duvIduv;
MisfitOuts.uAdjoint=[];
MisfitOuts.vAdjoint=[];



%% Calculate the gradient of the misfit function I with respect to the control variables (model parameters) p (here A and C or B).
%
% This is a bit tricky because I=I(u(p))
%
if CtrlVar.Inverse.CalcGradI
    
    
    switch lower(CtrlVar.Inverse.DataMisfit.GradientCalculation)
        
        case {"fixpoint","fixpointc","-fixpoint-","-fixpointc-"}
            
            switch CtrlVar.Inverse.InvertForField
                
                case "C"
                    
                    
                    DCI=FixPointGradHessianC(UserVar,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,RunInfo);
                    
                    
                case "B"
                    
                    
                    dBFuvLambda=Calc_FixPoint_deltaB(CtrlVar,MUA,F,Meas);
                    np=numel(dIdp); ddIdpp=sparse(np,np);
                    dBJ=0;
                    DBI=dBFuvLambda+dBJ;
                    
                otherwise
                    
                    fprintf(" CtrlVar.Inverse.InvertFor has an invalid value.\ n ")
                    fprintf(" CtrlVar.Inverse.InvertFor = %s \n ",CtrlVar.Inverse.InvertFor)
                    fprintf(" Fixpoint inversion only possible for C and B inversion. \n")
                    error("Misfit:IncorrectInputParameterCombination","Fixpoint inversion only possible for C and B inversion")
                    
            end
            
            
        case {"adjoint","-adjoint-"}
            %% Inverse problem
            %
            % Forward model:
            %   f(u(p),p)=0
            %
            %% Step 1: solve the non-linear forward problem (this has already be done ahead of the call to this m-file)
            %
            %       dfuv du = f(u)  ; u-> u+du until norm(f)<tolerance
            %
            %       [UserVar,RunInfo,F,l,drdu,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
            %
            %
            % [UserVar,RunInfo,F,l,dFduv,Ruv,Lubvb]= uv(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l);
            %% Step 2:  Solve adjoint equation, i.e.   dfuv l = -dJduv
            % fprintf(' Solve ajoint problem \n ')
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
            
            duvJ=duvIduv;     % Because this is the only J term that depends on (u,v). 
                              % If the regularization term depended on the state variable q, ie R=R(u,v) then this would not be correct.

            % Now solve the linear adjoint problem for lambda
            [lambda,lAdjoint]=solveKApeSymmetric(dfuv,LAdjoint,duvJ,LAdjointrhs,[],lAdjoint,CtrlVar);
            
            
            if CtrlVar.TestAdjointFiniteDifferenceType=="complex step differentiation"
                CtrlVar.TestForRealValues=false;
            end
            
            if CtrlVar.TestForRealValues && ~isreal(lAdjoint)
                save TestSave ; error("When solving adjoint equation Lagrange parmeters complex ")
            end
            
            uAdjoint=real(lambda(1:MUA.Nnodes)) ;
            vAdjoint=real(lambda(MUA.Nnodes+1:2*MUA.Nnodes));
            
            MisfitOuts.uAdjoint=uAdjoint;
            MisfitOuts.vAdjoint=vAdjoint;
            
            if CtrlVar.Inverse.InfoLevel>=1000 && CtrlVar.doplots
                
                GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,F.GF,CtrlVar);
                tri=MUA.connectivity;
                
                figure
                hold off
                subplot(2,2,1)
                [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,duvIduv(1:length(F.ub)),CtrlVar);  title("dIdu")
                hold on ; plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
                
                subplot(2,2,2)
                [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,duvIduv(1+length(F.ub):end),CtrlVar);  title('dIdv')
                hold on ; plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
                
                subplot(2,2,3)
                [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,uAdjoint,CtrlVar);  title('lx')
                hold on ; plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
                
                subplot(2,2,4)
                [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,vAdjoint,CtrlVar);  title('ly')
                hold on ; plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
            end
            
            %% Step 3:  <d_p F^* \lambda>,
            %
            % Note that I'm adding the d_p R term in the regularization step.
            %
            % But I need to include a possible <d_p I , \phi> term here.
            %
            % For p=A and p=C, d_p I =0 because I is not an explicit function of A and C
            %
            % But for b, d_b I = p_x (u db)
            
            if contains(lower(CtrlVar.Inverse.InvertFor),"c")
                
                dCFuvLambda=dIdCq(CtrlVar,UserVar,MUA,F,uAdjoint,vAdjoint,Meas);
                
                dCI=0 ;               % This is the explicit derivative of I with respect to C. 
                                      % The misfit term I is not an explicit function of C, so this equals to zero.

                DCI=dCFuvLambda+dCI;  % this is the part of the dI/dC derivative which is due to the implicit dependency 
                                      % of I on C because the velocities depend on C,

            end
            
            if contains(lower(CtrlVar.Inverse.InvertFor),"aglen")
                
               
                dAFuvLambda=dIdAq(CtrlVar,UserVar,MUA,F,uAdjoint,vAdjoint,Meas);
                
                dAI=0 ; % No explicit dependency of the misfit term I on A.

                DAI=dAFuvLambda+dAI;
            end
            
            
            if contains(CtrlVar.Inverse.InvertFor,"-B-")
                
                %  p= B ;
                
                dBdp=  1+zeros(MUA.Nnodes,1);
                %dBdp=  F.GF.node ; %
                dbdp=  F.GF.node ; % - (1-F.GF.node).*F.GF.node.*F.rho/F.rhow;
                dhdp= -F.GF.node ;
                
                % dIdB= dhF^* \lambda + dhJ
                % if only -dhdt- meas and no regularization
                % then dJdB=dh/db*dhJhdot
                
                dBFuvLambda=dIdbq(CtrlVar,MUA,uAdjoint,vAdjoint,F,dhdp,dbdp,dBdp);
             %   dBFuvLambda2=dIdBq2(CtrlVar,MUA,uAdjoint,vAdjoint,F);
                %dBFuvLambda=dBFuvLambda2;
                
                
                dBI=dhdp.*dhIhdot;  % The Ihdot misfit term includes an explicit dependency on B, which is here accounted for.
                DBI=dBFuvLambda+dBI;
                
                if CtrlVar.Inverse.OnlyModifyBedUpstreamOfGL
                    F.GF=IceSheetIceShelves(CtrlVar,MUA,F.GF,GLgeo,GLnodes,GLele) ;
                    DBI(~F.GF.NodesUpstreamOfGroundingLines)=0;
                end
                
                
                
            end
            
            
            
        otherwise
            
            fprintf(" CtrlVar.Inverse.DataMisfit.GradientCalculation has the value %s \n",CtrlVar.Inverse.DataMisfit.GradientCalculation)
            fprintf(" but the only allowed values are ''fixpoint'' and ''adjoint'' \n")
            error(" which case? ")
            
    end
 
    
    % Hessians
    
    if isfield(CtrlVar.Inverse.DataMisfit,'HessianEstimate')
        error(' field no longer used ')
    end
    
    if contains(CtrlVar.Inverse.MinimisationMethod,"Hessian")
        
        if contains(CtrlVar.Inverse.InvertForField,"C")
            
            if contains(CtrlVar.Inverse.Hessian,"IHC=FP")
                [~,ddIdCC]=FixPointGradHessianC(UserVar,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,RunInfo);
            elseif contains(CtrlVar.Inverse.Hessian,"IHC=GN")
                [ddIdCC]=GaussNewtonHessianC(UserVar,CtrlVar,MUA,DCI,F,Meas);
            elseif contains(CtrlVar.Inverse.Hessian,"IHC=M")
                ddIdCC=MUA.M/MUA.Area;
            elseif contains(CtrlVar.Inverse.Hessian,"IHC=D")
                ddIdCC=(MUA.Dxx+MUA.Dyy)/MUA.Area;
            elseif contains(CtrlVar.Inverse.Hessian,"IHC=0") || contains(CtrlVar.Inverse.Hessian,"IHC=O")
                N=MUA.Nnodes;
                ddIdCC=sparse(N,N);
            elseif  contains(CtrlVar.Inverse.Hessian,"IHC=I") || contains(CtrlVar.Inverse.Hessian,"IHC=1")
                N=MUA.Nnodes;
                ddIdCC=speye(N,N);
            else
                error('case not found')
            end
        end
        
        if contains(CtrlVar.Inverse.InvertForField,"A")
            if contains(CtrlVar.Inverse.Hessian,"IHA=FP")
                [~,ddIdAA]=FixPointGradHessianA(UserVar,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,RunInfo);
            elseif contains(CtrlVar.Inverse.Hessian,"IHA=GN")
                ddIdAA=GaussNewtonHessianA(UserVar,CtrlVar,MUA,DAI,F,Meas);
            elseif contains(CtrlVar.Inverse.Hessian,"IHA=M")
                ddIdAA=MUA.M/MUA.Area;
            elseif contains(CtrlVar.Inverse.Hessian,"IHA=D")
                ddIdAA=(MUA.Dxx+MUA.Dyy)/MUA.Area;
            elseif contains(CtrlVar.Inverse.Hessian,"IHA=0") || contains(CtrlVar.Inverse.Hessian,"IHA=O")
                N=MUA.Nnodes;
                ddIdAA=sparse(N,N);
            elseif contains(CtrlVar.Inverse.Hessian,"IHA=I") || contains(CtrlVar.Inverse.Hessian,"IHA=1")
                N=MUA.Nnodes;
                ddIdAA=speye(N,N);
            else
                error('case not found')
            end
        end
    end
    
    switch CtrlVar.Inverse.InvertForField
        
        case "A"
            dIdp=DAI;
            ddIdpp=ddIdAA ;
        case "b"
            error("fdsa")
        case "B"
            dIdp=DBI;
        case "C"
            dIdp=DCI;
            ddIdpp=ddIdCC ;
        case "AC"
            dIdp=[DAI;DCI];
            
            if contains(CtrlVar.Inverse.MinimisationMethod,"Hessian")
                N=MUA.Nnodes;
                ddIdpp = spalloc(N+N,N+N,nnz(ddIdAA)+nnz(ddIdCC));
                ddIdpp(1:N,1:N) = ddIdAA;
                ddIdpp(N+1:N+N,N+1:N+N) = ddIdCC;
            end
            
        otherwise
            
            error("sdfsa")
            
    end
    
    
else
    
    dIdp=0;
    
end

I=CtrlVar.Inverse.DataMisfit.Multiplier*I;

if nargout>1
    
    dIdp=CtrlVar.Inverse.DataMisfit.Multiplier*dIdp;
    ddIdpp=CtrlVar.Inverse.DataMisfit.Multiplier*ddIdpp;
    
end


if nargout>3
    MisfitOuts.I=I;
    MisfitOuts.dIdC=CtrlVar.Inverse.DataMisfit.Multiplier*DCI;
    MisfitOuts.dIdAGlen=CtrlVar.Inverse.DataMisfit.Multiplier*DAI;
    MisfitOuts.dIdB=CtrlVar.Inverse.DataMisfit.Multiplier*DBI;
end

end

