function Method=dhdtMisfitMethod(CtrlVar)

%%
%
%   Method=dhdtMisfitMethod(CtrlVar)
%
% Returns the formulation of the dh/dt misfit term, CtrlVar.Inverse.DataMisfit.dhdt.Method, as a string:
%
%   "-nodal-integral-" : (default) measurements at the nodes, misfit integrated over the mesh (EvaluateJhdotAndDerivatives.m)
%   "-points-"         : modelled dh/dt evaluated at the measurement locations, no averaging
%   "-blocks-"         : modelled dh/dt averaged over blocks (footprints) of the measurements
%
% If the field is not defined (for example CtrlVar from an older restart file) the default "-nodal-integral-" is returned.
%
% see also: dhdtObservationMisfit.m, CreateDhdtPointMeasurements.m, CreateDhdtBlockMeasurements.m
%
%%

Method="-nodal-integral-";

if isfield(CtrlVar,'Inverse') && isfield(CtrlVar.Inverse,'DataMisfit') && isfield(CtrlVar.Inverse.DataMisfit,'dhdt') ...
        && isfield(CtrlVar.Inverse.DataMisfit.dhdt,'Method')

    Method=string(CtrlVar.Inverse.DataMisfit.dhdt.Method);

    if ~ismember(Method,["-nodal-integral-","-points-","-blocks-"])
        error('dhdtMisfitMethod:UnknownMethod','CtrlVar.Inverse.DataMisfit.dhdt.Method=%s. Allowed values are -nodal-integral-, -points- and -blocks-.',Method)
    end

end

end
