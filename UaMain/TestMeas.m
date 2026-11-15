function Meas=TestMeas(CtrlVar,MUA,Meas)


[is_uv_meas,is_dhdt_meas]=is_uv_dhdt_Meas(CtrlVar);


if is_uv_meas && is_uv_points_misfit(CtrlVar)

    % "-points-" velocity misfit: the point measurements are used, and Meas.us, Meas.vs, Meas.usCov and Meas.vsCov are not (14 November 2026)
    if isempty(Meas.uvO) || isempty(Meas.usData) || isempty(Meas.vsData)
        fprintf(">>>>>>>>>>>>>>  CtrlVar.Inverse.DataMisfit.uv.Method=%s, but the point measurements of velocities are not defined!\n",CtrlVar.Inverse.DataMisfit.uv.Method)
        fprintf(" Define Meas.uvO, Meas.usData, Meas.vsData, Meas.usErr and Meas.vsErr using CreatePointMeasurements.m \n")
        error("TestMeas:MeasNotDefined","Point measurements of velocities not defined on input")
    end
    m=numel(Meas.usData);
    if ~isequal([numel(Meas.vsData) numel(Meas.usErr) numel(Meas.vsErr) size(Meas.uvO,1)],[m m m m])
        error("TestMeas:InconsistentPointMeasurements","Meas.usData, Meas.vsData, Meas.usErr, Meas.vsErr and Meas.uvO are not consistent in size.")
    end
    if size(Meas.uvO,2)~=MUA.Nnodes
        error("TestMeas:InconsistentPointMeasurements","Meas.uvO has %i columns, but the mesh has %i nodes. Has the mesh changed? Call CreatePointMeasurements.m again.",size(Meas.uvO,2),MUA.Nnodes)
    end
    if any(~isfinite([Meas.usData;Meas.vsData;Meas.usErr;Meas.vsErr])) || any([Meas.usErr;Meas.vsErr]<=0)
        error("TestMeas:InvalidPointMeasurements","Point measurements of velocities, or their errors, are not finite, or the errors are not strictly positive.")
    end

elseif is_uv_meas
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
        if any(diag(Meas.vsCov)==0)
            error('Ua:TestMeas','The error covariance matrix vsCov is a diagonal matrix and has zeros on the diagonal. This is not allowed. Modify data errors')
        end
    end
    
end


if is_dhdt_meas && dhdtMisfitMethod(CtrlVar)~="-nodal-integral-"

    % "-points-" or "-blocks-" dh/dt misfit: Meas.dhdt and Meas.dhdtCov are not used (14 November 2026)
    Method=dhdtMisfitMethod(CtrlVar);
    if isempty(Meas.dhdtData) || isempty(Meas.dhdtErr)
        fprintf(">>>>>>>>>>>>>>  CtrlVar.Inverse.DataMisfit.dhdt.Method=%s, but the dh/dt measurements are not defined!\n",Method)
        fprintf(" Define them using CreateDhdtPointMeasurements.m or CreateDhdtBlockMeasurements.m \n")
        error("TestMeas:MeasNotDefined","dh/dt measurements not defined on input")
    end
    m=numel(Meas.dhdtData);
    if numel(Meas.dhdtErr)~=m
        error("TestMeas:InconsistentDhdtMeasurements","Meas.dhdtData and Meas.dhdtErr are not consistent in size.")
    end
    if any(~isfinite([Meas.dhdtData(:);Meas.dhdtErr(:)])) || any(Meas.dhdtErr<=0)
        error("TestMeas:InvalidDhdtMeasurements","dh/dt measurements, or their errors, are not finite, or the errors are not strictly positive.")
    end
    switch Method
        case "-points-"
            if isempty(Meas.dhdtP) || numel(Meas.dhdtP.EleID)~=m || size(Meas.dhdtP.fun,2)~=MUA.nod || max(Meas.dhdtP.EleID)>MUA.Nele
                error("TestMeas:InconsistentDhdtMeasurements","Meas.dhdtP is not consistent with the measurements or the mesh. Has the mesh changed? Call CreateDhdtPointMeasurements.m again.")
            end
        case "-blocks-"
            if ~isequal(size(Meas.dhdtA),[m MUA.Nele*MUA.nip])
                error("TestMeas:InconsistentDhdtMeasurements","Meas.dhdtA is not consistent with the measurements, the mesh or the integration rule. Has the mesh changed? Call CreateDhdtBlockMeasurements.m again.")
            end
    end

elseif is_dhdt_meas
    
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
        
        if isscalar(Meas.dhdt)
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