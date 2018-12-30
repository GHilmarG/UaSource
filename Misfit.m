function [I,dIdp,ddIddp,MisfitOuts]=Misfit(UserVar,CtrlVar,MUA,BCs,F,l,GF,Priors,Meas,BCsAdjoint,RunInfo,dfuv)

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
%       2 I=d M d + d (Dx +Dy) d
%
%  where d=d_measured-d_calculated
%
%
Area=TriAreaTotalFE(MUA.coordinates,MUA.connectivity);

dIdC=[];
dIdAGlen=[];
dIdb=[];
dIdp=[] ;
ddIddp=sparse(1,1);

us=F.ub+F.ud ;
vs=F.vb+F.vd ;


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

if contains(CtrlVar.Inverse.Measurements,'-dhdt-','IgnoreCase',true)
        
        %%
        % CtrlVar.Tracer.SUPG.Use=1; CtrlVar.Tracer.SUPG.tau='tau2';
        % CtrlVar.dt=1e-6;
        % [UserVar,hnew,lambda]=TracerConservationEquation(UserVar,CtrlVar,MUA,CtrlVar.dt,F.h,F.ub,F.vb,F.as+F.ab,F.ub,F.vb,F.as+F.ab,0,BCs);
        % [hnew,l]=SSS2dPrognostic(CtrlVar,MUA,BCs,l,F.h,F.ub,F.vb,F.dubdt,F.dvbdt,F.as+F.ab,F.as*0,F.ub,F.vb,F.as+F.ab,F.as*0,F.dubdt,F.dvbdt);
        % [F.h hnew c1 c1-hnew]
        % F.dhdt=(hnew-F.h)/CtrlVar.dt;
        
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
        
        
        [UserVar,F.dhdt]=dhdtExplicit(UserVar,CtrlVar,MUA,F);
        
        if isdiag(Meas.dhdtCov)
            dhdtErr=sqrt(spdiags(Meas.dhdtCov));
            dhdtres=(F.dhdt-Meas.dhdt)./dhdtErr;
        else
            error('Misfit:Cov','Data covariance matrices must be diagonal')
        end
        
end

%% Calculate misfit term I and its gradient with respect to the state variables u and v
% This is straigtforward as the misfit term is an explicit function of u
% and v.
switch lower(CtrlVar.Inverse.DataMisfit.FunctionEvaluation)
    
    case 'uvdiscrete'  % this is here for comparision and testing, do not use, it's wrong!
        
        N=numel(us);
        I=full(usres'*usres+vsres'*vsres)/2/N;
        dIdu=usres./uErr/N;
        dIdv=vsres./vErr/N;
        dIduv=[dIdu(:);dIdv(:)];
        
    case 'integral' % always evaluate the continuous approximation
        
        if ~isfield(MUA,'M')
            MUA.M=MassMatrix2D1dof(MUA);
        end
        
        switch CtrlVar.Inverse.Measurements
            
            case '-uv-'
                
                dIdu=(MUA.M*usres)./uErr/Area;
                dIdv=(MUA.M*vsres)./vErr/Area;
                
                I=full(usres'*MUA.M*usres+vsres'*MUA.M*vsres)/2/Area;
                dIduv=[dIdu(:);dIdv(:)];
                
            case '-dhdt-'
                
                %dIdhdt=(MUA.M*dhdtres)./dhdtErr/Area;
                I=full(dhdtres'*MUA.M*dhdtres)/2/Area;
                
                [UserVar,dIhduv]=dIhdotduv(UserVar,CtrlVar,MUA,F,dhdtres,dhdtErr);
                dIhdu=dIhduv(1:MUA.Nnodes)/Area;
                dIhdv=dIhduv(MUA.Nnodes+1:end)/Area;
                dIduv=[dIhdu(:);dIhdv(:)];
                
            case {'-uv-dhdt-','-dhdt-uv-'}
                
                dIdu=(MUA.M*usres)./uErr/Area;
                dIdv=(MUA.M*vsres)./vErr/Area;
                %dIdhdt=(MUA.M*dhdtres)./dhdtErr/Area;
                [UserVar,dIhduv]=dIhdotduv(UserVar,CtrlVar,MUA,F,dhdtres,dhdtErr);
                dIhdu=dIhduv(1:MUA.Nnodes)/Area;
                dIhdv=dIhduv(MUA.Nnodes+1:end)/Area;
                
                I=full(usres'*MUA.M*usres+vsres'*MUA.M*vsres+dhdtres'*MUA.M*dhdtres)/2/Area;
                dIduv=[dIdu(:)+dIhdu(:);dIdv(:)+dIhdv(:)];
            otherwise
                
                fprintf('The case %s for the variable CtrlVar.Inverse.Measurements not recognized.\n',CtrlVar.Inverse.Measurements)
                error(' Redefine CtrlVar.Inverse.Measurements ')
                
        end
        
        if ~isreal(dIduv)  && CtrlVar.TestForRealValues
            save TestSave ; error('MisfitFunction:dIduvNoReal','dIduv is not real! Possibly a problem with covariance of data.')
        end
        
    otherwise
        error(' what case? ' )
end

MisfitOuts.dIduv=dIduv;
MisfitOuts.uAdjoint=[];
MisfitOuts.vAdjoint=[];


%% Calculate the gradient of the misfit function I with respect to the control variables (model parameters) p (here A and C).
%
% This is a bit tricky because I=I(u(p))
%
if CtrlVar.Inverse.CalcGradI
    
    
    switch lower(CtrlVar.Inverse.DataMisfit.GradientCalculation)
        
        case {'fixpoint','fixpointc','-fixpoint-','-fixpointc-'}
            
            if contains(lower(CtrlVar.Inverse.InvertFor),'c')
                
                if contains(lower(CtrlVar.Inverse.InvertFor),'aglen')
                    
                    fprintf(' CtrlVar.Inverse.InvertFor has an invalid value.\ n ')
                    fprintf(' CtrlVar.Inverse.InvertFor = %s \n ',CtrlVar.Inverse.InvertFor)
                    fprintf(' Fixpoint evaluation of data misfit gradient only possibly when solving for slipperiness alone. \n')
                    error('Fixpoint evaluation of data misfit gradient not implemented for AGlen or logAGlen inversion. ')
                    
                else
                    
                    dIdC=Calc_FixPoint_deltaC(CtrlVar,MUA,F.C,F.m,GF,F.ub,F.vb,Meas.us,Meas.vs);
                    np=numel(dIdp); ddIddp=sparse(np,np);
                    
                end
                
            else
                
                error('Fixpoint evaluation of data misfit gradient not implemented for AGlen or logAGlen inversion')
                
                
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
            
            
            dJdu=dIduv(:); % because the regularization term does not depend on u or v
            
            MLC_Adjoint=BCs2MLC(MUA,BCsAdjoint);
            LAdjoint=MLC_Adjoint.ubvbL;
            LAdjointrhs=MLC_Adjoint.ubvbRhs;
            lAdjoint=zeros(numel(LAdjointrhs),1) ;
            
            dJduAdjoint=dJdu;
            
            [lambda,lAdjoint]=solveKApeSymmetric(dfuv,LAdjoint,dJduAdjoint,LAdjointrhs,[],lAdjoint,CtrlVar);
            
            
            if ~isreal(lAdjoint)
                save TestSave ; error('When solving adjoint equation Lagrange parmeters complex ')
            end
            
            uAdjoint=real(lambda(1:MUA.Nnodes)) ; 
            vAdjoint=real(lambda(MUA.Nnodes+1:2*MUA.Nnodes));
            
            MisfitOuts.uAdjoint=uAdjoint;
            MisfitOuts.vAdjoint=vAdjoint;
            
            if CtrlVar.Inverse.InfoLevel>=1000 && CtrlVar.doplots
                
                GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
                tri=MUA.connectivity;
                
                figure
                hold off
                subplot(2,2,1)
                [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,dIduv(1:length(F.ub)),CtrlVar);  title('dIdu')
                hold on ; plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
                
                subplot(2,2,2)
                [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,dIduv(1+length(F.ub):end),CtrlVar);  title('dIdv')
                hold on ; plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
                
                subplot(2,2,3)
                [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,uAdjoint,CtrlVar);  title('lx')
                hold on ; plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
                
                subplot(2,2,4)
                [FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(tri,MUA.coordinates,vAdjoint,CtrlVar);  title('ly')
                hold on ; plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'r','LineWidth',2)
            end
            
            
            
            if contains(lower(CtrlVar.Inverse.InvertFor),'c')
                
                switch lower(CtrlVar.Inverse.DataGradient.FunctionEvaluation)
                    
                    case 'discrete' % Direct gradient evaluated from nodal points.
                        
                        if CtrlVar.CisElementBased
                            M = Ele2Nodes(MUA.connectivity,MUA.Nnodes);
                            Cnode=M*C;
                        else
                            Cnode=C;
                        end
                        
                        dIdC = -(1/m)*GF.node.*(Cnode+CtrlVar.CAdjointZero).^(-1/m-1).*(sqrt(ub.*ub+vb.*vb+CtrlVar.SpeedZero^2)).^(1/m-1).*(u.*uAdjoint+v.*vAdjoint);
                        
                        if contains(lower(CtrlVar.Inverse.InvertFor),'-logc-')
                            dIdC=log(10)*Cnode.*dIdp;
                        end
                        
                        if CtrlVar.CisElementBased
                            dIdC=Nodes2EleMean(MUA.connectivity,dIdC);
                        end
                        
                        np=numel(dIdp);
                        ddIddp=sparse(np,np);
                        
                    case 'integral'
                        
                        if CtrlVar.CisElementBased
                            dIdC=dIdCqEleSteps(CtrlVar,MUA,uAdjoint,vAdjoint,F.s,F.b,F.h,F.S,F.B,F.ub,F.vb,F.ud,F.vd,F.AGlen,F.n,F.C,F.m,F.rho,F.rhow,F.alpha,F.g,GF);
                            np=numel(dIdp); ddIddp=sparse(ones(np,1),1:np,1:np);
                        else
                            dIdC=dIdCq(CtrlVar,MUA,uAdjoint,vAdjoint,F.s,F.b,F.h,F.S,F.B,F.ub,F.vb,F.ud,F.vd,F.AGlen,F.n,F.C,F.m,F.rho,F.rhow,F.alpha,F.g,GF);
                        end
                end
            end
            
            if contains(lower(CtrlVar.Inverse.InvertFor),'aglen')
                
                switch lower(CtrlVar.Inverse.DataGradient.FunctionEvaluation)
                    
                    case 'discrete' % Direct gradient evaluated from nodal points.
                        
                        fprintf(' CtrlVar.AdjointGradientEvaluation=''uvdiscrete'' not possible in a combination with AGlen inversion.\n')
                        error('AdjointGradientNR2d:DiscreteAdjointAGlen','Discrete case not implemented. Used integral evaluation instead.')
                        
                    case 'integral'
                        
                        if CtrlVar.AGlenisElementBased
                            
                            dIdAGlen=dIdAEleSteps(CtrlVar,MUA,uAdjoint,vAdjoint,F.s,F.b,F.h,F.S,F.B,F.ub,F.vb,F.ud,F.vd,F.AGlen,F.n,F.C,F.m,F.rho,F.rhow,F.alpha,F.g,GF);
                            np=numel(dIdp); ddIddp=sparse(ones(np,1),1:np,1:np);
                            
                        else
                            dIdAGlen=dIdAq(CtrlVar,MUA,uAdjoint,vAdjoint,F.s,F.b,F.h,F.S,F.B,F.ub,F.vb,F.ud,F.vd,F.AGlen,F.n,F.C,F.m,F.rho,F.rhow,F.alpha,F.g,GF);
                        end
                end
            end
            
            
            if contains(lower(CtrlVar.Inverse.InvertFor),'-b-')
                
                switch lower(CtrlVar.Inverse.DataGradient.FunctionEvaluation)
                    
                    case 'discrete' % Direct gradient evaluated from nodal points.
                        
                        fprintf(' CtrlVar.AdjointGradientEvaluation=''uvdiscrete'' not possible in a combination with b inversion.\n')
                        error('AdjointGradientNR2d:DiscreteAdjointAGlen','Discrete case not implemented. Used integral evaluation instead.')
                        
                    case 'integral'
                        
                        dIdb=dIdbq(CtrlVar,MUA,uAdjoint,vAdjoint,F);
                        
                end
            end
            
            
            
        otherwise
            
            fprintf(' CtrlVar.Inverse.DataMisfit.GradientCalculation has the value %s \n',CtrlVar.Inverse.DataMisfit.GradientCalculation)
            fprintf(' but the only allowed values are ''fixpoint'' and ''adjoint'' \n')
            error(' which case? ')
            
    end
    
    
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
            dIdp=dIdAGlen;
        case 'b'
            dIdp=dIdb;
        case 'C'
            dIdp=dIdC;
        case 'Ab'
            dIdp=[dIdAGlen;dIdb];
        case 'AC'
            dIdp=[dIdAGlen;dIdC];
        case 'bC'
            dIdp=[dIdb;dIdC];
        case 'AbC'
            dIdp=[dIdAGlen;dIdb;dIdC];
            
        otherwise
            
            dIdp=0;
            error('sdfsa')
            
    end
    
    %     if contains(lower(CtrlVar.Inverse.InvertFor),'aglen') && contains(lower(CtrlVar.Inverse.InvertFor),'c')
    %
    %         dIdp=[dIdAGlen;dIdC];
    %
    %     elseif contains(lower(CtrlVar.Inverse.InvertFor),'aglen') && ~contains(lower(CtrlVar.Inverse.InvertFor),'c')
    %
    %         dIdp=dIdAGlen;
    %
    %     elseif ~contains(lower(CtrlVar.Inverse.InvertFor),'aglen') && contains(lower(CtrlVar.Inverse.InvertFor),'c')
    %
    %         dIdp=dIdC;
    %
    %     elseif ~contains(lower(CtrlVar.Inverse.InvertFor),'aglen') && contains(lower(CtrlVar.Inverse.InvertFor),'c')
    %
    %     end
    
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
    MisfitOuts.dIdC=CtrlVar.Inverse.DataMisfit.Multiplier*dIdC;
    MisfitOuts.dIdAGlen=CtrlVar.Inverse.DataMisfit.Multiplier*dIdAGlen;
    MisfitOuts.dIdb=CtrlVar.Inverse.DataMisfit.Multiplier*dIdb;
end




end

