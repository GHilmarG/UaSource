function [J,dJduv,ObjFunc,RunInfo,UserVar]=MisfitFunction(UserVar,RunInfo,CtrlVar,MUA,F,l,Priors,Meas)

%%
% J is the object function
% ObjFunc is a structure with fields giving all terms of the object function,
% i.e misfit term, regularistation terms, penalty terms etc. 
% 
% J=ObjFunc.I+ObjFunc.R
%
% In the optimisation J and dJduv are used, but not ObjFunc directly. 
%
%
%%



%[J,Idata,IRegC,IRegAGlen,dIduv,IBarrierC,IBarrierAGlen]=MisfitFunction(UserVar,CtrlVar,MUA,ub,vb,ud,vd,AGlen,C,Priors,Meas)

narginchk(11,11)

us=F.ub+F.ud ; vs=F.vb+F.vd; % Surface velocities

%     J=Idata+IRegC+IRegAGlen+IBarrierC+IBarrierAGlen;
%
%%
%
% The misfit function is 
%
% $$J=I + R + B$$ 
%
% where $I$ is a data misfit term, $R$ a regularisation term, and $B$ a barrier
% term.
%
% The data mistfit term is on the form
%
% $$ I=\frac{1}{2\cal{A}} \int \int (u(x,y)-\hat{u}(x,y)) \, \kappa(x-x',y-y') \, (u(x',y')-\hat{u}(x',y')) dx \, dy \, dx' \, dy' $$
%
% where $\kappa$ is a (inverse) covariance kernel, $u$ are model output quantity, and $\hat{u}$ measurements of that quantity. 
%
% It is assumed that data errors are spatially uncorrelated, i.e.
%
% $$\kappa(x,y)=e^{-2} \; \delta(x,y)$$ 
%
% where $e$ are data errors. Therefore
%
% $$ I=\frac{1}{2\cal{A}} \int (u(x,y)-\hat{u}(x,y))/e)^2 dx \, dy $$
%
% $I$ is dimentionless.
%
% In FE context:
%
% $$ I= \frac{1}{2\cal{A}}  ((\mathbf{u}-\hat{\mathbf{u}})./\mathbf{e})' \, \mathbf{M} \,((\mathbf{u}-\hat{\mathbf{u}})./\mathbf{e}) $$
%
% here $\mathbf{u}$ etc, are vectors and we used the ./ matlab notation to
% indicate element wise division.
% 
% The derivative of the data mistfit term with respect ot $u$ is
%
% $$ dI/du =  \frac{1}{\cal{A}}  \mathbf{M} \,(\mathbf{u}-\hat{\mathbf{u}}) $$
%
%
% The regularisation term $R$ is written in descrete form as
%
% $$R= \frac{1}{2 N} (\mathbf{c}-\hat{\mathbf{c}})' \mathbf{K}^{-1} (\mathbf{c}-\hat{\mathbf{c}})$$
%
% where $\mathbf{}$ is a covariance matrix and $N$ the number of elements in $\mathbf{c}$. $R$ is dimentionless
%
% $$dR/dc=\frac{1}{N} \mathbf{K}^{-1} (\mathbf{c}-\hat{\mathbf{c}})$$
%
%%
Area=TriAreaTotalFE(MUA.coordinates,MUA.connectivity);

% Calculate data misfit functions and gradients with respect to observables u
% and v.
switch lower(CtrlVar.MisfitFunction)
    
    case 'uvdiscrete'
        
        N=numel(us);
        Idata=(us-Meas.us)'*(Meas.usCov\(us-Meas.us))/(2*N)+(vs-Meas.vs)'*(Meas.vsCov\(vs-Meas.vs))/(2*N);
        dJdu=(1/N)*Meas.usCov\(us-Meas.us) ;
        dJdv=(1/N)*Meas.vsCov\(vs-Meas.vs);
        dJduv=[dJdu(:) ; dJdv(:)];
        
    case 'uvintegral'

         if isdiag(Meas.usCov)
             uErr=sqrt(spdiags(Meas.usCov)); 
             usres=(us-Meas.us)./uErr;
         else
             error('MisfitFunction:Cov','Data covariance matrices but be diagonal')
         end
         
         if isdiag(Meas.vsCov)
             vErr=sqrt(spdiags(Meas.vsCov)); 
             vsres=(vs-Meas.vs)./vErr;
         else
             error('MisfitFunction:Cov','Data covariance matrices must be diagonal')
         end
        
        M=MassMatrix2D1dof(MUA);
        Idata=full(usres'*M*usres/2+vsres'*M*vsres/2)/Area;
        dJdu=(M*usres)./uErr/Area;
        dJdv=(M*vsres)./vErr/Area;
        dJduv=[dJdu(:);dJdv(:)];
        
        if ~isreal(dJduv)
            save TestSave ; error('MisfitFunction:dJduvNoReal','dJduv is not real! Possibly a problem with covariance of data.')
        end
        
                
    otherwise
        error(' what case? ' )
end

ObjFunc.I=CtrlVar.MisfitMultiplier*Idata;


ObjFunc.RegC=Calc_IRegC(CtrlVar,MUA,Priors.CovC,F.C,Priors.C);
ObjFunc.RegAGlen=Calc_IRegdAGlen(CtrlVar,MUA,Priors.CovAGlen,F.AGlen,Priors.AGlen);
ObjFunc.BarrierC=Calc_IBarrierC(CtrlVar,F.C);
ObjFunc.BarrierAGlen=Calc_IBarrierAGlen(CtrlVar,F.AGlen);

ObjFunc.R=ObjFunc.RegC+ObjFunc.RegAGlen+ObjFunc.BarrierC+ObjFunc.BarrierAGlen;

J=ObjFunc.I+ObjFunc.R;  


%fprintf('MisfitFunction: J=%-g \t Idata=%-g \t IRegC=%-g \t IRegAGlen=%-g \t IBarrierC=%-g \t IBarrierAGlen=%-g \n',...
%    J,Idata,IRegC,IRegAGlen,IBarrierC,IBarrierAGlen)
end

