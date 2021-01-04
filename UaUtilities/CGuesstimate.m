function [UserVar,CGuess,DrivingStress,GradSurf]=CGuesstimate(UserVar,CtrlVar,MUA,F,GF,Meas,CorrelationDistance)

%%
%
%   [UserVar,CGuess,DrivingStress]=CGuesstimate(UserVar,CtrlVar,MUA,F,GF,Meas,CorrelationDistance)
%
% Provides a rough estimate of C based on measured surface velocities and
% an estimate of driving stress assuming Weertman sliding law. Can easily be
% modifed to include Budd or other Robin type BCs.
%
% Returns
%
%   C=speed/DrivingStress^m
%
% Optionally, smooths estimate using Helmholtz smoothing for a given Matern correlation
% distance.
%
% Example:
%
%   CorrelationDistance=50e3; 
%   [UserVar,CGuess,DrivingStress,GradSurf]=CGuesstimate(UserVar,CtrlVar,MUA,F,F.GF,Meas,CorrelationDistance); 
%
%%


DrivingStressMax=200;
DrivingStressMin=5;

GradSurfMin=1e-4; 
GradSurfMax=0.1; 
BasalSpeedMax=3000 ;
CGuessMax=mean(BasalSpeedMax./(DrivingStressMin.^F.m));

% Calc surface gradient
[dsdx,dsdy]=calcFEderivativesMUA(F.s,MUA,CtrlVar);
[dsdx,dsdy]=ProjectFintOntoNodes(MUA,dsdx,dsdy);



GradSurf=sqrt(dsdx.*dsdx+dsdy.*dsdy) ;
GradSurf(GradSurf<GradSurfMin)=GradSurfMin ;
GradSurf(GradSurf>GradSurfMax)=GradSurfMax ;

SpeedMeasured=sqrt(Meas.us.*Meas.us+Meas.vs.*Meas.vs) ;
DrivingStress= F.rho.*F.g.*(F.h+CtrlVar.ThickMin).*GradSurf;

DrivingStress(DrivingStress>DrivingStressMax)=DrivingStressMax;
DrivingStress(DrivingStress<DrivingStressMin)=DrivingStressMin;

CGuess=SpeedMeasured./(DrivingStress.^F.m);
CGuess(CGuess>CGuessMax)=CGuessMax;

%  u = c rho g h ds/dx


if nargin>6 && ~isnan(CorrelationDistance)
    % Helmholtz smoothing
    % This is really just here as an example
    % Arguably better to do whatever smoothing required outside of this routine
    Dimention=2; alpha=2 ; % get the wavenumber in the Helmoltz equation corresponding to a given
    % correlation distance for the Matern covariance of degree one
    
    
    [~,~,kappa,~]=Matern(alpha,CorrelationDistance,Dimention,NaN);
    
    [UserVar,GradSurf]=HelmholtzEquation(UserVar,CtrlVar,MUA,1,1/kappa^2,GradSurf,0);
    
    GradSurf(GradSurf<GradSurfMin)=GradSurfMin ;
    GradSurf(GradSurf>GradSurfMax)=GradSurfMax ;
    
    
    DrivingStress= F.rho.*F.g.*(F.h+CtrlVar.ThickMin).*GradSurf;
    DrivingStress(DrivingStress>DrivingStressMax)=DrivingStressMax;
    DrivingStress(DrivingStress<DrivingStressMin)=DrivingStressMin;
    
    % Smoothing driving stress
    [UserVar,DrivingStress]=HelmholtzEquation(UserVar,CtrlVar,MUA,1,1/kappa^2,DrivingStress,0);
    DrivingStress(DrivingStress>DrivingStressMax)=DrivingStressMax;
    DrivingStress(DrivingStress<DrivingStressMin)=DrivingStressMin;
    CGuess=SpeedMeasured./(DrivingStress.^F.m);
    CGuess(CGuess>CGuessMax)=CGuessMax;
end


end



