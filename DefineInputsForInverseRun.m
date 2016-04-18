function [InvStartValues,Priors,Meas,BCsAdjoint]=DefineInputsForInverseRun(Experiment,CtrlVar,MUA,BCs,InvStartValues,Priors,Meas,BCsAdjoint,time,AGlen,C,n,m,s,b,S,B,rho,rhow,GF)

%%%
% An example of how to define inputs for an inverse run
%
% If doing an inverse run then a good starting point might be to 
% put a copy of this file into your local run directory and then change as needed.
%
%

x=MUA.coordinates(:,1) ; y=MUA.coordinates(:,2);

%% define measurements and measurement errors
% Here it is assumed that one is inverting measurements of surface velocities 
% and that scattered interpolants have already been defined and stored in a file
% 
fprintf(' Loading measured velocities \n')
load('MyMeasurements','Fu','Fv')

Meas.us=Fu(x,y);  % Mapping measurments onto FE mesh
Meas.vs=Fv(x,y);
Meas.ws=Meas.vs*0 ;

% Defining data erros and data covariance matrices.
% Often errors are uncorelated and the covariance matrices
% are therefore diagonal.
usError=1; vsError=1  ; wsError=1;
Meas.usCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,usError.^2,MUA.Nnodes,MUA.Nnodes);
Meas.vsCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,vsError.^2,MUA.Nnodes,MUA.Nnodes);
Meas.wsCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,wsError.^2,MUA.Nnodes,MUA.Nnodes);

% if add errors
%uMeas=uMeas+uError.*randn(MUA.Nnodes,1);
%vMeas=vMeas+vError.*randn(MUA.Nnodes,1);
%wMeas=wMeas+wError*randn(MUA.Nnodes,1);

%% define boundary conditions of adjoint problem
% The boundary conditions of the adjoint problem need to be defined.
% Often the BCs of the adjoint problem are homogeneous.
% If the forward problem uses periodic boundary conditions,
% then those of the adjoint problem are usually periodic as well

% If, for example, adjoint problem has the same boundary conditions as the forward, set: 
% BCsAdjoint=BCs; % this would, for example, be the case if the forward problem has periodic BCs,
                  % in which case the adjoint problem will have periodic BCs as well.

% If the adjoint problem has homogenous BCs (a very common case) then do: 
BCsAdjoint.ubFixedNode=MUA.Boundary.Nodes ; BCsAdjoint.ubFixedValue=BCsAdjoint.ubFixedNode*0;
BCsAdjoint.vbFixedNode=MUA.Boundary.Nodes ; BCsAdjoint.vbFixedValue=BCsAdjoint.vbFixedNode*0;



%% now define priors
% covariance matrices for priors
% these covariance matrices are typically not diagonal
% make sure to define those on nodes if A and C are defined on nodes
% and on elements if A and C are defined as elemetn values.
% This here is just an examle and might need to be adjusted.
if CtrlVar.AGlenisElementBased
    CAGlen=sparse(1:MUA.Nele,1:MUA.Nele,1,MUA.Nele,MUA.Nele);
else
    CAGlen=sparse(1:MUA.Nnodes,1:MUA.Nnodes,1,MUA.Nnodes,MUA.Nnodes);
end

if CtrlVar.isRegC    % only create a full covariance matrix if isRegC is on.
    Err=1e-5 ; Sigma=10e3 ; DistanceCutoff=10*Sigma;
    
    if CtrlVar.CisElementBased
        xEle=mean(reshape(x(MUA.connectivity,1),MUA.Nele,MUA.nod),2);  yEle=mean(reshape(y(MUA.connectivity,1),MUA.Nele,MUA.nod),2);
        [CC]=SparseCovarianceDistanceMatrix(xEle,yEle,Err,Sigma,DistanceCutoff);
    else
        [CC]=SparseCovarianceDistanceMatrix(x,y,Err,Sigma,DistanceCutoff);
    end
    
else
    if CtrlVar.CisElementBased
        CC=sparse(1:MUA.Nele,1:MUA.Nele,1,MUA.Nele,MUA.Nele);
    else
        CC=sparse(1:MUA.Nnodes,1:MUA.Nnodes,1,MUA.Nnodes,MUA.Nnodes);
    end
end

Priors.CovAGlen=CAGlen;
Priors.CovC=CC;

Priors.s=s;
Priors.b=b;
Priors.S=S;
Priors.B=B;

% Priors for C and A need to be defined
% Here something simple is prescribed but this might be based on some
% previous estimates, for example A could be based on some temperature model, etc.
% 
Priors.C=1e-6+zeros(MUA.Nnodes,1); % assuming C defined on nodes
Priors.m=3;
Priors.AGlen=AGlenVersusTemp(-10)+zeros(MUA.Nnodes,1); 
Priors.n=3;

Priors.rho=rho;
Priors.rhow=rhow;

%% Define start values
% I'm here setting starting values equal to priors
[InvStartValues.C,InvStartValues.m]=DefineSlipperyDistribution(Experiment,CtrlVar,MUA,time,s,b,s-b,S,B,rho,rhow,GF);
[InvStartValues.AGlen,InvStartValues.n]=DefineAGlenDistribution(Experiment,CtrlVar,MUA,time,s,b,s-b,S,B,rho,rhow,GF);


end
