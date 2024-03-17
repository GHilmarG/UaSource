function [UserVar,F]=GetCalving(UserVar,CtrlVar,MUA,F,BCs)


narginchk(5,5)
nargoutchk(2,2)


if ~CtrlVar.LevelSetMethod
    return
end

nArgs=nargin('DefineCalving');



switch nArgs

    case 5

        [UserVar,F.LSF,F.c]=DefineCalving(UserVar,CtrlVar,MUA,F,BCs) ;

    case 7

        [UserVar,F.LSF,F.c]=DefineCalving(UserVar,CtrlVar,MUA,F.LSF,F.c,F,BCs) ;


    otherwise

        error('DefineCalving must have either 5 or 7 inputs arguments.')

end

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
%     if any(isnan(F.c))  && CtrlVar.CurrentRunStepNumber>0
%         errorStruct.identifier = 'GetCalving:NaNinInput';
%         errorStruct.message = 'nan in c (calving)';
%         error(errorStruct)
%     end


if numel(F.c)==1
    F.c=F.c+zeros(MUA.Nnodes,1);
end

if CtrlVar.LevelSetEvolution=="-prescribed-"
    F.c=nan;
end


if numel(F.LSF)==1
    F.LSF=F.LSF+zeros(MUA.Nnodes,1);
end

F.LSFMask=CalcMeshMask(CtrlVar,MUA,F.LSF,0);



% else
%     if ~isempty(F.c)  &&  numel(F.c)~=MUA.Nnodes
%         errorStruct.identifier = 'GetCalving:CalvingFieldInvalid';
%         errorStruct.message = 'number of elements in the calving field must equal number of nodes';
%         error(errorStruct)
%     end
% end
%% Modify c away from calving front
%
%
% The calving rate can be effectivly nautralized by setting c=(u,v) grad \varphi
% The automated mass-balance feedback is applied to LSFMask.NOdesOut so I must




end