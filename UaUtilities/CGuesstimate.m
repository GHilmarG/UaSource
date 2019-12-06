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
%%

% Calc surface gradient

[dsdx,dsdy]=calcFEderivativesMUA(F.s,MUA,CtrlVar);
[dsdx,dsdy]=ProjectFintOntoNodes(MUA,dsdx,dsdy);

if ~isfield(CtrlVar,'SurfaceSlopeMin')
    CtrlVar.SurfaceSlopeMin=1e-4;
end


GradSurf=sqrt(dsdx.*dsdx+dsdy.*dsdy) ;
Speed=sqrt(Meas.us.*Meas.us+Meas.vs.*Meas.vs) ;
DrivingStress= F.rho.*F.g.*(F.h+CtrlVar.ThickMin).*(GradSurf+CtrlVar.SurfaceSlopeMin);

CGuess=Speed./(DrivingStress.^F.m);


%  u = c rho g h ds/dx


if nargin>8 && ~isnan(CorrelationDistance)
    % Helmholtz smoothing
    % This is really just here as an example
    % Arguably better to do whatever smoothing required outside of this routine
    Dimention=2; alpha=2 ; % get the wavenumber in the Helmoltz equation corresponding to a given
                           % correlation distance for the Matern covariance of degree one
                           
        
    [~,~,kappa,~]=Matern(alpha,CorrelationDistance,Dimention,NaN);
    
    % Smoothing driving stress
    [UserVar,DrivingStress]=HelmholtzEquation(UserVar,CtrlVar,MUA,1,1/kappa^2,DrivingStress,0);
    
    % Smoothing surface speed
    [UserVar,Speed]=HelmholtzEquation(UserVar,CtrlVar,MUA,1,1/kappa^2,Speed,0);
    CGuess=Speed./(DrivingStress.^F.m);
    
    % Can't stop smoothing...
    [UserVar,CGuess]=HelmholtzEquation(UserVar,CtrlVar,MUA,1,1/kappa^2,CGuess,0);
end


end



