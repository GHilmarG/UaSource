function Meas=TestMeas(CtrlVar,MUA,Meas)


if any(isnan(Meas.us)) ; error('Ua:TestMeas:NaN','NaN in surface measurements of u') ; end

if any(isnan(Meas.vs)) ; error('Ua:TestMeas:NaN','NaN in surface measurements of v') ; end

if any(isnan(Meas.ws)) ; error('Ua:TestMeas:NaN','NaN in surface measurements of w') ; end


if isdiag(Meas.usCov)
    if any(diag(Meas.usCov)==0)
        error('Ua:TestMeas','The error covariance matrix usCov is a diagonal matrix and has zeros on the diagonal. This is not allowed. Modify data errors')
    end
end


if isdiag(Meas.vsCov)
    if any(diag(Meas.usCov)==0)
        error('Ua:TestMeas','The error covariance matrix usCov is a diagonal matrix and has zeros on the diagonal. This is not allowed. Modify data errors')
    end
end

if isdiag(Meas.wsCov)
    if any(diag(Meas.usCov)==0)
        error('Ua:TestMeas','The error covariance matrix usCov is a diagonal matrix and has zeros on the diagonal. This is not allowed. Modify data errors')
    end
end

end