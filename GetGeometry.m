function [s,b,S,B,alpha]=GetGeometry(Experiment,CtrlVar,MUA,time,FieldsToBeDefined)


[s,b,S,B,alpha]=DefineGeometry(Experiment,CtrlVar,MUA,time,FieldsToBeDefined);


% some error checks
errorStruct.identifier = 'GetGeometry:NaNinInput';

if ~isempty(strfind(FieldsToBeDefined,'s'))
    if any(isnan(s))
        errorStruct.message = 'nan in s';
        error(errorStruct)
    end
end

if ~isempty(strfind(FieldsToBeDefined,'b'))
    if any(isnan(b))
        errorStruct.message = 'nan in b';
        error(errorStruct)
    end
end

if ~isempty(strfind(FieldsToBeDefined,'S'))
    if any(isnan(S))
        errorStruct.message = 'nan in S';
        error(errorStruct)
    end
end

if ~isempty(strfind(FieldsToBeDefined,'B'))
    if any(isnan(B))
        errorStruct.message = 'nan in B';
        error(errorStruct)
    end
end




end