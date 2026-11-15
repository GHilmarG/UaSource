function [Meas,Info]=CreatePointMeasurements(CtrlVar,MUA,Meas,x,y,us,vs,usErr,vsErr)

%%
%
%   [Meas,Info]=CreatePointMeasurements(CtrlVar,MUA,Meas,x,y,us,vs,usErr,vsErr)
%
% Defines the point measurements of surface velocities used by the "-points-" velocity misfit,
%
%   CtrlVar.Inverse.DataMisfit.uv.Method="-points-" ;
%
% in which the modelled velocities are evaluated at the measurement locations (the mesh is brought to the data), in the same way as
% for the direct measurements of B. The misfit is
%
% $$ I_{uv} = \frac{1}{2} \sum_k \frac{((\mathbf{O}\mathbf{u})_k-u_k)^2}{\sigma_{u,k}^2} + \frac{1}{2} \sum_k \frac{((\mathbf{O}\mathbf{v})_k-v_k)^2}{\sigma_{v,k}^2} $$
%
% (see PointDataMisfit.m), where O is the node-to-data operator (see BuildNode2DataMap.m).
%
% Inputs:
%
%   x, y          : locations of the measurements, for example the pixel centres of a gridded velocity product
%   us, vs        : measured velocity components at (x,y)
%   usErr, vsErr  : measurement errors (standard deviations) of us and vs
%
% Measurements outside the mesh, and measurements where any of x, y, us, vs, usErr, vsErr is not finite, or where an error is not
% strictly positive, are discarded. On return, the following fields of Meas are defined (one entry per retained measurement):
%
%   Meas.uvx, Meas.uvy        : locations
%   Meas.usData, Meas.vsData  : measured values
%   Meas.usErr, Meas.vsErr    : errors (standard deviations)
%   Meas.uvO                  : node-to-data operator, numel(Meas.uvx) x MUA.Nnodes
%
% Note: Meas.uvO depends on the mesh. If the mesh changes, this function must be called again.
%
% *Errors*
%
% The misfit treats the errors as uncorrelated. If the measurement errors are spatially correlated, this can be accounted for in an
% approximate way when defining the errors here, for example by thinning the data (ThinMeasurementsOnGrid.m), or by using an
% effective variance. For a gridded product with pixel area a_p, uncorrelated (nugget) error sigma_p, and spatially correlated errors
% with variance sigma_M^2 and a Matern correlation (nu=1) with range rho, the effective variance is
%
% $$ \sigma_{\mathrm{eff}}^2 = \sigma_p^2 + \sigma_M^2 \, \frac{A_{\mathrm{eff}}}{a_p} , \qquad A_{\mathrm{eff}}=\frac{\pi}{2} \rho^2 $$
%
% This is appropriate if the elements are larger than rho, and a_p < A_eff. Here a_p is the area per (retained) measurement.
%
% *Example*
%
%   [Meas,Info]=CreatePointMeasurements(CtrlVar,MUA,Meas,xPix,yPix,uPix,vPix,uErrPix,vErrPix) ;
%
% see also: PointDataMisfit.m, BuildNode2DataMap.m, is_uv_points_misfit.m
%
%%

narginchk(9,9)

x=x(:) ; y=y(:) ; us=us(:) ; vs=vs(:) ; usErr=usErr(:) ; vsErr=vsErr(:) ;

n=numel(x);
if ~all([numel(y) numel(us) numel(vs) numel(usErr) numel(vsErr)]==n)
    error("CreatePointMeasurements:InconsistentInputs","x, y, us, vs, usErr and vsErr must all have the same number of elements.")
end

Valid=isfinite(x) & isfinite(y) & isfinite(us) & isfinite(vs) & isfinite(usErr) & isfinite(vsErr) & usErr>0 & vsErr>0 ;
iValid=find(Valid);

[O,Inside]=BuildNode2DataMap(CtrlVar,MUA,x(iValid),y(iValid));

iKeep=iValid(Inside);

Meas.uvx=x(iKeep);
Meas.uvy=y(iKeep);
Meas.usData=us(iKeep);
Meas.vsData=vs(iKeep);
Meas.usErr=usErr(iKeep);
Meas.vsErr=vsErr(iKeep);
Meas.uvO=O(Inside,:);

Info.nInput=n;
Info.nNotValid=n-numel(iValid);
Info.nOutsideMesh=numel(iValid)-numel(iKeep);
Info.nRetained=numel(iKeep);

if isempty(CtrlVar) || ~isfield(CtrlVar,"InfoLevel") || CtrlVar.InfoLevel>=1
    fprintf("CreatePointMeasurements: %i velocity measurements retained (of %i; %i not valid, %i outside the mesh).\n", ...
        Info.nRetained,Info.nInput,Info.nNotValid,Info.nOutsideMesh)
end

if Info.nRetained==0
    warning("CreatePointMeasurements:NoData","No velocity measurements retained.")
end

end
