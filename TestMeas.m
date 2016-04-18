function Meas=TestMeas(CtrlVar,MUA,Meas)


if any(isnan(Meas.us)) ; error('Ua:TestMeas:NaN','NaN in surface measurements of u') ; end

if any(isnan(Meas.vs)) ; error('Ua:TestMeas:NaN','NaN in surface measurements of v') ; end

if any(isnan(Meas.ws)) ; error('Ua:TestMeas:NaN','NaN in surface measurements of w') ; end


end