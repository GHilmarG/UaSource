function [I,dIdp,ddIddp,MisfitOuts]=Misfit(UserVar,CtrlVar,MUA,BCs,F,l,Priors,Meas,BCsAdjoint,RunInfo,dfuv)

%%
%
%           J(u,p) = I(u(p)) + R(p)
%
% Forward model:   f(u(p),p)=0
%
%
% Calculate data misfit functions and gradients with respect to control
% variables p and the state variables u and v.
%
%   dfuv    : The derivative of the forward model f(u(p),p)=0 (scalar) with respect to
%             the state variables (uv).
%
%
%   As explained in Ua compendium the misfit term, I, can be written as
%
%       2 I=d M d + d (Dx +Dy) d%
%  where d=d_measured-d_calculated
%
%
%  Notation:
%  J = Juv   + Jdot    + JA    +   JB   + JC
%
%  I need:
%
%       duvJ  =   duvJuv   + duvJdot    (i.e. any J terms that are function of uv)
%       dCJC
%       dBJB
%       dAJA
%
%
%   Adjoint Eqs:  duvFuv lambda = duvJ
%
%       DAJ     = dAFuv lambda + dAJ
%       DBJ     = dBFuv lambda + dBJ
%       DCJ     = dCFuv lambda + dCJ
%


persistent GLgeo GLnodes GLele

Area=TriAreaTotalFE(MUA.coordinates,MUA.connectivity);

dAFuvLambda=[];
dBFuvLambda=[];
dCFuvLambda=[];

DAJ=[];
DBJ=[];
DCJ=[];

dIdp=[] ;
ddIddp=sparse(1,1);

MisfitOuts.dIduv=[];
MisfitOuts.uAdjoint=[];
MisfitOuts.vAdjoint=[];



us=F.ub+F.ud ;
vs=F.vb+F.vd ;

%% Do some test on inputs, check if error covariance matrices are correctly defined
% calculate residual terms, i.e. difference between meas and calc values.
if contains(CtrlVar.Inverse.Measurements,'-uv-','IgnoreCase',true)
    
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

%% Data misfit terms and derivatives
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

%% Calculate misfit term I and its gradient with respect to the state variables u and v
% This is straigtforward as the misfit term is an explicit function of u
% and v.

Juv=0 ; Jhdot=0 ;
duJdu=zeros(MUA.Nnodes,1) ;
dvJdv=zeros(MUA.Nnodes,1) ;

duJhdot=zeros(MUA.Nnodes,1);
dvJhdot=zeros(MUA.Nnodes,1);
dhJhdot=zeros(MUA.Nnodes,1) ;

switch lower(CtrlVar.Inverse.DataMisfit.FunctionEvaluation)
    
    case 'uvdiscrete'  % this is here for comparision and testing, do not use, it's wrong!
        
        N=numel(us);
        Juv=full(usres'*usres+vsres'*vsres)/2/N;
        duJdu=usres./uErr/N;
        dvJdv=vsres./vErr/N;
        
        
    case 'integral' % always evaluate the continuous approximation
        
        if ~isfield(MUA,'M')
            MUA.M=MassMatrix2D1dof(MUA);
        end
        
        if contains(CtrlVar.Inverse.Measurements,'-uv-')
            
            duJdu=(MUA.M*usres)./uErr/Area;
            dvJdv=(MUA.M*vsres)./vErr/Area;
            Juv=full(usres'*MUA.M*usres+vsres'*MUA.M*vsres)/2/Area;
            
        end
        
        
        if contains(CtrlVar.Inverse.Measurements,'-dhdt-')
            
            [Jhdot,duJhdot,dvJhdot,dhJhdot]=EvaluateJhdotAndDerivatives(UserVar,CtrlVar,MUA,F,BCs,Meas);
            
        end
        
        
        
end

I=Juv+Jhdot ;  %  but still missing the regularisation terms: JA, JB and JC
duvJduv=[duJdu(:)+duJhdot(:);dvJdv(:)+dvJhdot(:)];

if CtrlVar.Inverse.TestAdjoint.FiniteDifferenceType=="complex step differentiation"
    CtrlVar.TestForRealValues=false;
end

if ~isreal(duvJduv)  && CtrlVar.TestForRealValues
    save TestSave ; error('MisfitFunction:dIduvNoReal','dIduv is not real! Possibly a problem with covariance of data.')
end


MisfitOuts.dIduv=duvJduv;
MisfitOuts.uAdjoint=[];
MisfitOuts.vAdjoint=[];



%% Calculate the gradient of the misfit function I with respect to the control variables (model parameters) p (here A and C or B).
%
% This is a bit tricky because I=I(u(p))
%
if CtrlVar.Inverse.CalcGradI
    
    
    switch lower(CtrlVar.Inverse.DataMisfit.GradientCalculation)
        
        case {'fixpoint','fixpointc','-fixpoint-','-fixpointc-'}
            
            switch CtrlVar.Inverse.InvertForField
                
                case 'C'
                    
                    %dCFuvLambda=Calc_FixPoint_deltaC(CtrlVar,MUA,F.C,F.m,F.GF,F.ub,F.vb,Meas.us,Meas.vs);
                    dCFuvLambda=Calc_FixPoint_deltaC(CtrlVar,UserVar,MUA,F,Meas);
                    np=numel(dIdp); ddIddp=sparse(np,np);
                    dCJ=0 ;
                    DCJ=dCFuvLambda+dCJ;
                    
                case 'B'
                    
                    
                    dBFuvLambda=Calc_FixPoint_deltaB(CtrlVar,MUA,F,Meas);
                    np=numel(dIdp); ddIddp=sparse(np,np);
                    dBJ=0;
                    DBJ=dBFuvLambda+dBJ;
                    
                otherwise
                    
                    fprintf(' CtrlVar.Inverse.InvertFor has an invalid value.\ n ')
                    fprintf(' CtrlVar.Inverse.InvertFor = %s \n ',CtrlVar.Inverse.InvertFor)
                    fprintf(' Fixpoint inversion only possible for C and B inversion. \n')
                    error('Misfit:IncorrectInputParameterCombination','Fixpoint inversion only possible for C and B inversion')
                    
            end
            
            
        case 'adjoint'
            %% Inverse problem
            %
            % Forward model:
            %   f(u(p),p)=0
            %
            %% Step 1: solve the non-linear forward problem
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
            
            duvJ=duvJduv;     % Because this is the only J terms that depents on UV
            [lambda,lAdjoint]=solveKApeSymmetric(dfuv,LAdjoint,duvJ,LAdjointrhs,[],lAdjoint,CtrlVar);
            
            
            if CtrlVar.Inverse.TestAdjoint.FiniteDifferenceType=="complex step differentiation"
                CtrlVar.TestForRealValues=false;
            end
            
            if CtrlVar.TestForRealValues && ~isreal(lAdjoint)
                save TestSave ; error('When solving adjoint equation Lagrange parmeters complex ')
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
                [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,duvJduv(1:length(F.ub)),CtrlVar);  title('dIdu')
                hold on ; plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
                
                subplot(2,2,2)
                [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,duvJduv(1+length(F.ub):end),CtrlVar);  title('dIdv')
                hold on ; plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
                
                subplot(2,2,3)
                [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,uAdjoint,CtrlVar);  title('lx')
                hold on ; plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
                
                subplot(2,2,4)
                [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,vAdjoint,CtrlVar);  title('ly')
                hold on ; plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
            end
            
            %% Step 3:  <d_p F^* \lambda>,
            % Note that I'm adding the d_p R term in the regularisation step
            % But I need to include a possible <d_p I , \phi> term here
            % For p=A and p=C, d_p I =0 because I is not an explicit function of A and C
            % But for b, d_b I = p_x (u db)
            
            if contains(lower(CtrlVar.Inverse.InvertFor),'c')
                
                switch lower(CtrlVar.Inverse.DataGradient.FunctionEvaluation)
                    
                    case 'discrete' % Direct gradient evaluated from nodal points.
                        
                        if CtrlVar.CisElementBased
                            M = Ele2Nodes(MUA.connectivity,MUA.Nnodes);
                            Cnode=M*F.C;
                        else
                            Cnode=F.C;
                        end
                        
                        dCFuvLambda = -(1./F.m)*F.GF.node.*(Cnode+CtrlVar.CAdjointZero).^(-1./F.m-1).*(sqrt(F.ub.*F.ub+F.vb.*F.vb+CtrlVar.SpeedZero^2)).^(1./F.m-1).*(F.ub.*uAdjoint+F.vb.*vAdjoint);
                        
                        if contains(lower(CtrlVar.Inverse.InvertFor),'-logc-')
                            dCFuvLambda=log(10)*Cnode.*dCFuvLambda;
                        end
                        
                        if CtrlVar.CisElementBased
                            dCFuvLambda=Nodes2EleMean(MUA.connectivity,dCFuvLambda);
                        end
                        
                        np=numel(dCFuvLambda);
                        ddIddp=sparse(np,np);
                        
                    case 'integral'
                        
                        if CtrlVar.CisElementBased
                            dCFuvLambda=dIdCqEleSteps(CtrlVar,MUA,uAdjoint,vAdjoint,F.s,F.b,F.h,F.S,F.B,F.ub,F.vb,F.ud,F.vd,F.AGlen,F.n,F.C,F.m,F.rho,F.rhow,F.alpha,F.g,F.GF);
                            np=numel(F.C); ddIddp=sparse(ones(np,1),1:np,1:np);
                        else
                            
                             dCFuvLambda=dIdCq(CtrlVar,UserVar,MUA,F,uAdjoint,vAdjoint);
                            %dCFuvLambda=dIdCq(CtrlVar,MUA,uAdjoint,vAdjoint,F.s,F.b,F.h,F.S,F.B,F.ub,F.vb,F.ud,F.vd,F.AGlen,F.n,F.C,F.m,F.rho,F.rhow,F.alpha,F.g,F.GF);
                        end
                end
                dCJ=0 ; %  Here I should add the regularisation term, rather then doing this outside of this function
                DCJ=dCFuvLambda+dCJ;
            end
            
            if contains(lower(CtrlVar.Inverse.InvertFor),'aglen')
                
                switch lower(CtrlVar.Inverse.DataGradient.FunctionEvaluation)
                    
                    case 'discrete' % Direct gradient evaluated from nodal points.
                        
                        fprintf(' CtrlVar.AdjointGradientEvaluation=''uvdiscrete'' not possible in a combination with AGlen inversion.\n')
                        error('AdjointGradientNR2d:DiscreteAdjointAGlen','Discrete case not implemented. Used integral evaluation instead.')
                        
                    case 'integral'
                        
                        if CtrlVar.AGlenisElementBased
                            
                            dAFuvLambda=dIdAEleSteps(CtrlVar,MUA,uAdjoint,vAdjoint,F.s,F.b,F.h,F.S,F.B,F.ub,F.vb,F.ud,F.vd,F.AGlen,F.n,F.C,F.m,F.rho,F.rhow,F.alpha,F.g,F.GF);
                            np=numel(dIdp); ddIddp=sparse(ones(np,1),1:np,1:np);
                            
                        else
                            dAFuvLambda=dIdAq(CtrlVar,MUA,uAdjoint,vAdjoint,F.s,F.b,F.h,F.S,F.B,F.ub,F.vb,F.ud,F.vd,F.AGlen,F.n,F.C,F.m,F.rho,F.rhow,F.alpha,F.g,F.GF);
                        end
                end
                dAJ=0 ; %  Here I should add the regularisation term, rather then doing this outside of this function
                DAJ=dAFuvLambda+dAJ;
            end
            
            
            if contains(CtrlVar.Inverse.InvertFor,'-B-')
                
                switch lower(CtrlVar.Inverse.DataGradient.FunctionEvaluation)
                    
                    case 'discrete' % Direct gradient evaluated from nodal points.
                        
                        fprintf(' CtrlVar.AdjointGradientEvaluation=''uvdiscrete'' not possible in a combination with b inversion.\n')
                        error('AdjointGradientNR2d:DiscreteAdjointAGlen','Discrete case not implemented. Used integral evaluation instead.')
                        
                    case 'integral'
                        
                        %  p= B ;
                        
                        dBdp=  1+zeros(MUA.Nnodes,1);
                        %dBdp=  F.GF.node ; %
                        dbdp=  F.GF.node ; % - (1-F.GF.node).*F.GF.node.*F.rho/F.rhow;
                        dhdp= -F.GF.node ;
                        
                        % dIdB= dhF^* \lambda + dhJ
                        % if only -dhdt- meas and no regularisation
                        % then dJdB=dh/db*dhJhdot
                        
                        dBFuvLambda=dIdbq(CtrlVar,MUA,uAdjoint,vAdjoint,F,dhdp,dbdp,dBdp);
                        dBFuvLambda2=dIdBq2(CtrlVar,MUA,uAdjoint,vAdjoint,F);
                        %dBFuvLambda=dBFuvLambda2; 
                        
                        
                        
                end
                dBJ=dhdp.*dhJhdot;
                DBJ=dBFuvLambda+dBJ;
                
                if CtrlVar.Inverse.OnlyModifyBedUpstreamOfGL
                    [F.GF,GLgeo,GLnodes,GLele]=IceSheetIceShelves(CtrlVar,MUA,F.GF,GLgeo,GLnodes,GLele) ;
                    DBJ(~F.GF.NodesUpstreamOfGroundingLines)=0;
                end
                

                
            end
            
            
            
        otherwise
            
            fprintf(' CtrlVar.Inverse.DataMisfit.GradientCalculation has the value %s \n',CtrlVar.Inverse.DataMisfit.GradientCalculation)
            fprintf(' but the only allowed values are ''fixpoint'' and ''adjoint'' \n')
            error(' which case? ')
            
    end
    
    
%figure ; PlotMeshScalarVariable(CtrlVar,MUA,DBJ) ; title('DBJ')
% [I,L,U,C] = isoutlier(DBJ,'gesd'); factor=1 ; DBJ(DBJ>(factor*U))=factor*U ; DBJ(DBJ<(L/factor))=L/factor ; 
%figure ; PlotMeshScalarVariable(CtrlVar,MUA,DBJ) ; title('DBJ')
    
% Hessians
    
    switch lower(CtrlVar.Inverse.InvertFor)
        
        case {'c','-c-'}
            Hscale=1/(mean(F.C)^2);
        case {'aglen','-aglen-'}
            Hscale=1/(mean(F.AGlen)^2);
        case {'logc','-logc-'}
            Hscale=1/(log10(mean(F.C)^2));
        case {'logaglen','-logaglen-'}
            Hscale=1/(log10(mean(F.AGlen)^2));
    end
    
    switch upper(CtrlVar.Inverse.DataMisfit.HessianEstimate)
        
        
        case {'0','O'}
            np=numel(dIdp);
            ddIddp=sparse(np,np);
        case {'1','I'}
            np=numel(dIdp);
            ddIddp=Hscale*sparse(ones(np,1),1:np,1:np);
        case 'M'
            ddIddp=Hscale*MUA.M;
    end
    
    
    switch CtrlVar.Inverse.InvertForField
        
        case 'A'
            dIdp=DAJ;
        case 'b'
            error('fdsa')
        case 'B'
            dIdp=DBJ;
        case 'C'
            dIdp=DCJ;
        case 'AC'
            dIdp=[DAJ;DCJ];
        otherwise
            
            error('sdfsa')
            
    end
    
    
else
    
    dIdp=0;
    
end

I=CtrlVar.Inverse.DataMisfit.Multiplier*I;

if nargout>1
    
    dIdp=CtrlVar.Inverse.DataMisfit.Multiplier*dIdp;
    ddIddp=CtrlVar.Inverse.DataMisfit.Multiplier*ddIddp;
    
end


if nargout>3
    MisfitOuts.I=I;
    MisfitOuts.dIdC=CtrlVar.Inverse.DataMisfit.Multiplier*DCJ;
    MisfitOuts.dIdAGlen=CtrlVar.Inverse.DataMisfit.Multiplier*DAJ;
    MisfitOuts.dIdB=CtrlVar.Inverse.DataMisfit.Multiplier*DBJ;
end

end

