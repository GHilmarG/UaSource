function [UserVar,BCsLevelSet,F]=GetCalving(UserVar,CtrlVar,MUA,BCs,BCsLevelSet,F)
    
    
    narginchk(6,6)
    nargoutchk(3,3)
    
    
    
    [UserVar,F.LSF,BCsLevelSet,F.c]=DefineCalving(UserVar,CtrlVar,MUA,BCs,BCsLevelSet,F) ;
    
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