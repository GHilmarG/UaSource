
function [is_uv_meas,is_dhdt_meas]=is_uv_dhdt_Meas(CtrlVar)


if contains(CtrlVar.Inverse.Measurements,"-uv-","IgnoreCase",true)
    is_uv_meas=true;
else
    is_uv_meas=false;
end



if contains(CtrlVar.Inverse.Measurements,'-dhdt-','IgnoreCase',true)
    is_dhdt_meas=true;
else
    is_dhdt_meas=false;
end


end