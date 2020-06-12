function [UserVar,F]=GetCalving(UserVar,CtrlVar,MUA,F,BCs)
    
    
    narginchk(5,5)
    nargoutchk(2,2)
    
    
    if ~CtrlVar.LevelSetMethod
        return
    end
    
    if CtrlVar.CalvingLaw=="-No Ice Shelves-"
        if isempty(F.LSF) && ~isempty(F.GF.node)
            F.LSF=ReinitializeLevelSet([],[],CtrlVar,MUA,F.GF.node,CtrlVar.GLthreshold);
        end
        return
    end
    
     
    [UserVar,F.LSF,F.c]=DefineCalving(UserVar,CtrlVar,MUA,F,BCs) ;
    
    % some input checks
    
    
    if any(isnan(F.LSF))
        errorStruct.identifier = 'GetCalving:NaNinInput';
        errorStruct.message = 'nan in LSF';
        error(errorStruct)
    end
    
    % I do allow LSF not being updated in each call, this may imply that LSF has wrong number
    % of elements, but this will be picked up in the mapping from old to new mesh
    %
    %     if numel(F.LSF)~=MUA.Nnodes
    %         errorStruct.identifier = 'GetCalving:LSFinvalid';
    %         errorStruct.message = 'number of elements in LSF must equal number of nodes';
    %         error(errorStruct)
    %     end
    %
    
%     if numel(F.LSF)~=MUA.Nnodes
%         fprintf(' Just checking take this out \n ' ) 
%     end
%     
    if any(isnan(F.c))
        errorStruct.identifier = 'GetCalving:NaNinInput';
        errorStruct.message = 'nan in c (calving)';
        error(errorStruct)
    end
    
    if numel(F.c)==1
        F.c=F.c+zeros(MUA.Nnodes,1);
    end
    
    if numel(F.c)~=MUA.Nnodes
        errorStruct.identifier = 'GetCalving:CalvingFieldInvalid';
        errorStruct.message = 'number of elements in the calving field must equal number of nodes';
        error(errorStruct)
    end
    
end