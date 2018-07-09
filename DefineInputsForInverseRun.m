function [UserVar,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo]=...
    DefineInputsForInverseRun(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,Priors,Meas,BCsAdjoint,RunInfo)

%%
% *Note: This m-file is just an example of how to define inputs for an inverse run. You will need to modify to fit your own problem.*
%
% What you need to define are:
% 
% 
% # Measurments and data errors (data errors are specified as diagonal covariance matrices.)
% # Start values for inversion. (These are some values for the model parameters that you want to invert for.)
% # Priors for the inverted fields. (Currently the only priors that are used the the priors for C and AGlen.)
%
%
% Note: When doing an inverse run, presumably a good start is to copy this file from the source directory to you own run director. 
%

persistent FuMeas FvMeas FerrMeas  % keep scattered interpolants for the data in memory. 


%% Load measurments and define data errors as diagonal covariance matrices
if isempty(FuMeas)  % Have I already loaded the data, if so don't do it again.
    
    fprintf('Loading interpolants for surface velocity data: %-s ',UserVar.SurfaceVelocityInterpolant)
    load(UserVar.SurfaceVelocityInterpolant,'FuMeas','FvMeas','FerrMeas')
    fprintf(' done.\n')
end

Meas.us=double(FuMeas(MUA.coordinates(:,1),MUA.coordinates(:,2)));
Meas.vs=double(FvMeas(MUA.coordinates(:,1),MUA.coordinates(:,2)));
Err=double(FerrMeas(MUA.coordinates(:,1),MUA.coordinates(:,2)));

MissingData=isnan(Meas.us) | isnan(Meas.vs) | isnan(Err) | (Err==0);
Meas.us(MissingData)=0 ;  Meas.vs(MissingData)=0 ; Err(MissingData)=1e10; 

usError=Err ; vsError=Err ; wsError=usError*0+1e10;
Meas.usCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,usError.^2,MUA.Nnodes,MUA.Nnodes);
Meas.vsCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,vsError.^2,MUA.Nnodes,MUA.Nnodes);
Meas.wsCov=sparse(1:MUA.Nnodes,1:MUA.Nnodes,wsError.^2,MUA.Nnodes,MUA.Nnodes);

%% Define start values for inversion. These are some values for the model parameters that you want to invert for

% Here the values for AGlen and C are obtained by calling these m-files. Typically, however the start values would be some reasonable initial
% estimates, or the results of a previous inversion. 
[UserVar,InvStartValues.C,InvStartValues.m]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.s-F.b,F.S,F.B,F.rho,F.rhow,GF);
[UserVar,InvStartValues.AGlen,InvStartValues.n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,CtrlVar.time,F.s,F.b,F.s-F.b,F.S,F.B,F.rho,F.rhow,GF);




%% Define Priors. Currently the only priors that are used the the priors for C and AGlen.
Priors.s=F.s;
Priors.b=F.b;
Priors.S=F.S;
Priors.B=F.B;
Priors.AGlen=AGlenVersusTemp(-10);
Priors.n=3; 
Priors.m=3; ub=10 ; tau=80 ; % units meters, year , kPa
C0=ub/tau^Priors.m;
Priors.C=C0;
Priors.rho=F.rho;
Priors.rhow=F.rhow;




%% Covariance of prior (if using Bayesian Regularisation)
% listingCC=dir('CC.mat') ; listingCA=dir('CAGlen.mat') ;
%
% 
%  Note: this is only used if using Bayesian Regularisation involving covariance matrices and when doing an element-wise inversion.
%
%  By default inversion is done over nodes and using Tikhonov regularisation. Hence, this defining covariance matrices for the priors is
%  then not needed. 
%
% if strcmpi(CtrlVar.Inverse.Regularize.Field,'cov')
%     CreateCovMatAndSave=1;
%     if numel(listingCC)==1 && numel(listingCA)==1
%         CreateCovMatAndSave=0;
%         FileName='CC.mat';
%         fprintf('DefineInverseModelingParameters: loading CC from file: %-s ',FileName)
%         load(FileName,'CC') ;
%         fprintf(' done \n ')
%         %%
%         
%         FileName='CAGlen.mat';
%         fprintf('DefineInverseModelingParameters: loading CAGlen from file: %-s ',FileName)
%         load(FileName,'CAGlen');
%         fprintf(' done \n ')
%         
%         if length(CC)~=length(F.C)
%             CreateCovMatAndSave=1;
%             fprintf(' Covariance matrix in input file does not have correct dimentions. Will create a new one \n')
%         end
%     end
%     
%     if CreateCovMatAndSave
%         
%         
%         xEle=Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,1)); yEle=Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,2));
%         Err=1e-1 ; Sigma=200 ; DistanceCutoff=5*Sigma;
%         fprintf('creating sparse covariance matrix ')  ; tStart=tic;
%         [CC]=SparseCovarianceDistanceMatrix(xEle,yEle,Err,Sigma,DistanceCutoff);
%         tElapsed=toc(tStart);
%         fprintf('in %-g sec \n',tElapsed)
%         FileName='CC.mat'; save(FileName,'CC')
%         
%         Err=1e-5 ; Sigma=200 ; DistanceCutoff=5*Sigma;
%         fprintf('creating sparse covariance matrix ')  ; tStart=tic;
%         [CAGlen]=SparseCovarianceDistanceMatrix(xEle,yEle,Err,Sigma,DistanceCutoff);
%         tElapsed=toc(tStart);
%         fprintf('in %-g sec \n',tElapsed)
%         FileName='CAGlen.mat'; save(FileName,'CAGlen')
%     end
% else
%     CC=[] ;
%     CAGlen=[];
% end
% 
% Priors.CovAGlen=CAGlen;
% Priors.CovC=CC;
% 
% 



end
