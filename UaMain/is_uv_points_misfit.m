function isPoints=is_uv_points_misfit(CtrlVar)

%%
%
%   isPoints=is_uv_points_misfit(CtrlVar)
%
% Returns true if the velocity misfit is evaluated using point measurements, i.e. with the modelled velocities evaluated at the
% measurement locations (CtrlVar.Inverse.DataMisfit.uv.Method="-points-"), and false if the (default) nodal-integral formulation is
% used (CtrlVar.Inverse.DataMisfit.uv.Method="-nodal-integral-").
%
% If the field is not defined (for example CtrlVar from an older restart file) the default "-nodal-integral-" is assumed.
%
% see also: PointDataMisfit.m, CreatePointMeasurements.m, Ua2D_DefaultParameters.m
%
%%

isPoints=false;

if isfield(CtrlVar,"Inverse") && isfield(CtrlVar.Inverse,"DataMisfit") && isfield(CtrlVar.Inverse.DataMisfit,"uv") ...
        && isfield(CtrlVar.Inverse.DataMisfit.uv,"Method")

    switch string(CtrlVar.Inverse.DataMisfit.uv.Method)
        case "-points-"
            isPoints=true;
        case "-nodal-integral-"
            isPoints=false;
        otherwise
            error("is_uv_points_misfit:UnknownMethod","CtrlVar.Inverse.DataMisfit.uv.Method=%s. Allowed values are -nodal-integral- and -points-.",string(CtrlVar.Inverse.DataMisfit.uv.Method))
    end

end

end
