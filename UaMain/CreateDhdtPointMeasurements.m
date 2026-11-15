function [Meas,Info]=CreateDhdtPointMeasurements(CtrlVar,MUA,Meas,x,y,dhdt,dhdtErr)

%%
%
%   [Meas,Info]=CreateDhdtPointMeasurements(CtrlVar,MUA,Meas,x,y,dhdt,dhdtErr)
%
% Defines measurements of dh/dt at the locations (x,y) for the "-points-" dh/dt misfit,
%
%   CtrlVar.Inverse.DataMisfit.dhdt.Method="-points-" ;
%
% in which the modelled dh/dt is evaluated at the measurement locations, without any averaging (see dhdtObservationMisfit.m):
%
% $$ J = \frac{1}{2} \sum_k \frac{(\dot{h}(\mathbf{x}_k) - \dot{h}_k)^2}{\sigma_k^2} $$
%
% Measurements outside the mesh, not finite, or with errors that are not strictly positive, are discarded. Any further selection, for
% example only measurements over grounded ice, must be done before calling this function.
%
% On return, Meas.dhdtx, Meas.dhdty, Meas.dhdtData, Meas.dhdtErr and Meas.dhdtP (form functions and their derivatives at the
% measurement locations, see ShapeFunctionsAtPoints.m) are defined. Meas.dhdtP depends on the mesh, and this function must be called
% again if the mesh changes.
%
% Errors are treated as uncorrelated (see the comments on effective variances in CreatePointMeasurements.m).
%
% see also: CreateDhdtBlockMeasurements.m, dhdtObservationMisfit.m, ShapeFunctionsAtPoints.m
%
%%

narginchk(7,7)

x=x(:) ; y=y(:) ; dhdt=dhdt(:) ; dhdtErr=dhdtErr(:) ;
n=numel(x);
if ~all([numel(y) numel(dhdt) numel(dhdtErr)]==n)
    error('CreateDhdtPointMeasurements:InconsistentInputs','x, y, dhdt and dhdtErr must all have the same number of elements.')
end

Valid=isfinite(x) & isfinite(y) & isfinite(dhdt) & isfinite(dhdtErr) & dhdtErr>0 ;
iValid=find(Valid);

P=ShapeFunctionsAtPoints(CtrlVar,MUA,x(iValid),y(iValid));
iKeep=iValid(P.Inside);

Meas.dhdtx=x(iKeep);
Meas.dhdty=y(iKeep);
Meas.dhdtData=dhdt(iKeep);
Meas.dhdtErr=dhdtErr(iKeep);
Meas.dhdtP=rmfield(P,'Inside');
Meas.dhdtA=[];
Meas.dhdtBlockSize=[];

Info.nInput=n;
Info.nNotValid=n-numel(iValid);
Info.nOutsideMesh=numel(iValid)-numel(iKeep);
Info.nRetained=numel(iKeep);

if isempty(CtrlVar) || ~isfield(CtrlVar,'InfoLevel') || CtrlVar.InfoLevel>=1
    fprintf('CreateDhdtPointMeasurements: %i dh/dt measurements retained (of %i; %i not valid, %i outside the mesh).\n', ...
        Info.nRetained,Info.nInput,Info.nNotValid,Info.nOutsideMesh)
end

if Info.nRetained==0
    warning('CreateDhdtPointMeasurements:NoData','No dh/dt measurements retained.')
end

end
