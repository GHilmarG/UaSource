function OutsideValue=DefineOutsideValues(UserVar,CtrlVar,MUA,F,OutsideValue)

%%
%
%   OutsideValue=DefineOutsideValues(UserVar,CtrlVar,MUA,F,OutsideValue)
%
% Defines values outside of a previous computational domain.
%
% These values are used in a transient run when mapping
% (interpolating/extrapolating) from one mesh to another whenever some nodes of
% the new mesh are outside of the previous mesh. That is, when the new mesh has
% some nodes that are not witin any of the elements of the previous mesh.
%
% This situation can arise when using, for example, manual deactivation of
% elements. If the domain increases in size, then some intial velocites and
% thicknesses must be defined over these new areas.
%
% Typically, some minimum thickness will be defined over the new areas and the
% velocites set to zero. The velocity values are only used as a start values for
% a diagnostic uv solution. Hence, the exact velocites prescibed are not
% important (provided the uv solution converges). However, the thickness
% prescibed over the new areas is important as it defines the new start values
% over any outside areas. Note that in a diagnostic (time-independent run) the
% thickness is always prescribed using DefineGeometry.m.  For that reason this
% routine is never used in a diagnostic run.
%
% Only return scalar values.
%
%%

% Prescribe here outside thickness values of twice the minimum ice thickness.
OutsideValue.h=2*CtrlVar.ThickMin ;

% Make sure the s and b correspondes to flotation. However this is not essential
% as s and b are always adjusted internally based on h, S and B given rho and
% rhow.
OutsideValue.s=mean(F.S)+OutsideValue.h*(1-mean(F.rho)/F.rhow);
OutsideValue.b=OutsideValue.s-OutsideValue.h;

% Define reasonably initial values for velocities for the uv solver. These
% values will have no impact on the final solution (provided the non-linear uv
% solver converges using these initial estimates.).
OutsideValue.ub=0;
OutsideValue.vb=0;

OutsideValue.ud=0;
OutsideValue.vd=0;

% these rates are only needed for explicit estimates of velocities. Similarly to
% the uv values above, these exact values will have no impact on the solution
% (again provided the non-linear solver converges using these initial
% estimates.).
OutsideValue.dubdt=0;
OutsideValue.dvbdt=0;

OutsideValue.duddt=0;
OutsideValue.dvddt=0;

OutsideValue.dhdt=0;


end
