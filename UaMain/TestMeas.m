function Meas=TestMeas(CtrlVar,MUA,Meas)


if contains(CtrlVar.Inverse.Measurements,"-uv-")
    if isempty(Meas.us)
        fprintf(">>>>>>>>>>>>>>  Meas.us is empty!\n")
        fprintf(" CtrlVar.Inverse.Measurements=%s \n",CtrlVar.Inverse.Measurements)
        fprintf(" Meas.us needs to be defined. \n")
        error("TestMeas:MeasNotDefined","Measurments not defined on input")
    end
    
    if any(isnan(Meas.us)) ; error('Ua:TestMeas:NaN','NaN in surface measurements of u') ; end
    if any(isnan(Meas.vs)) ; error('Ua:TestMeas:NaN','NaN in surface measurements of v') ; end
    
    
    if isdiag(Meas.usCov)
        if any(diag(Meas.usCov)==0)
            error('Ua:TestMeas','The error covariance matrix usCov is a diagonal matrix and has zeros on the diagonal. This is not allowed. Modify data errors')
        end
    end
    
    
    if isdiag(Meas.vsCov)
        if any(diag(Meas.usCov)==0)
            error('Ua:TestMeas','The error covariance matrix vsCov is a diagonal matrix and has zeros on the diagonal. This is not allowed. Modify data errors')
        end
    end
    
end


if contains(CtrlVar.Inverse.Measurements,"-dhdt-")
    
    if isempty(Meas.dhdt)
        fprintf(">>>>>>>>>>>>>>  Meas.dhdt is empty!\n")
        fprintf(" CtrlVar.Inverse.Measurements=%s \n",CtrlVar.Inverse.Measurements)
        fprintf(" Meas.dhdt needs to be defined. \n")
        error("TestMeas:MeasNotDefined","Measurments not defined on input")
    end
    if any(isnan(Meas.dhdt)) ; error('Ua:TestMeas:NaN','NaN in surface measurements of dh/dt') ; end
    
    if isempty(Meas.dhdtCov)
        fprintf(">>>>>>>>>>>>>>  Meas.dhdtCov is empty!\n")
        fprintf(" CtrlVar.Inverse.Measurements=%s \n",CtrlVar.Inverse.Measurements)
        fprintf(" Meas.dhdtCov needs to be defined. \n")
        error("TestMeas:MeasNotDefined","Measurments not defined on input")
    end
    
    
    if isdiag(Meas.dhdtCov)
        if any(diag(Meas.dhdtCov)==0)
            error('Ua:TestMeas','The error covariance matrix dhdtCov is a diagonal matrix and has zeros on the diagonal. This is not allowed. Modify data errors')
        end
    end
    
    if numel(Meas.dhdt) ~= MUA.Nnodes
        
        if numel(Meas.dhdt)==1
            Meas.dhdt=Meas.dhdt+zeros(MUA.Nnodes,1);
        else
            fprintf(' Meas.dhdt must have same number of elements as the number of nodes.')
            fprintf(' numel(Meas.dhdt)=%i \n',numel(Meas.dhdt))
            fprintf(' MUA.Nnodes=%i \n',MUA.Nnodes)
            error(' Meas.dhdt has incorrect dimentions')
        end
        
    end
    
end


end