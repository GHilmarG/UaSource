





function [UserVar,RunInfo,F,l,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened]=GetDesiredEleSize(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened,NodalErrorIndicators)


persistent nCalls

narginchk(11,11)
nargoutchk(7,7)

InputFile="DefineDesiredEleSize.m"; TestIfInputFileInWorkingDirectory(InputFile) ;


if isempty(nCalls)
    nCalls=0;
end
nCalls=nCalls+1;


nOut=nargout('DefineDesiredEleSize');
nIn=nargin('DefineDesiredEleSize') ;

if nOut~=7 && nIn~=13

    fprintf("DefineDesiredEleSize.m does not have the right number of input or output arguments.\")

    fprintf(" The call to DefineDesiredEleSize must be on the form: \n ")
    fprintf("%s \n","[UserVar,RunInfo,F,l,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened]DefineDesiredEleSize(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,x,y,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened,NodalErrorIndicators)")

    error("Incorrect input.\ n")

end


if contains(lower(CtrlVar.MeshRefinementMethod),"local")
    x=MUA.xEle ; y=MUA.yEle;
else
    x=MUA.coordinates(:,1) ;
    y=MUA.coordinates(:,2) ;
end

InputEleSizeDesired=EleSizeDesired;
InputElementsToBeRefined=ElementsToBeRefined;


[UserVar,RunInfo,F,l,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened]=...
    DefineDesiredEleSize(UserVar,RunInfo,CtrlVar,MUA,BCs,F,l,x,y,EleSizeDesired,ElementsToBeRefined,ElementsToBeCoarsened,NodalErrorIndicators);



if contains(CtrlVar.MeshRefinementMethod,'local','IgnoreCase',true)

    if ~isequal(InputEleSizeDesired,EleSizeDesired)
        if nCalls==1
            fprintf(['Note: The mesh refinement methods is %s. \n' ...
                'When using this local mesh-refinement method in combination with ``DefineDesiredEleSize'' \n',...
                'you only need to specify as an output which elements are to be refined and/or coarsended. \n',...
                'You can specify ``EleSizeDesired'' as well, but it will be ignored.\n'],...
                CtrlVar.MeshRefinementMethod)
        end
    end

    % just some basic check here, not covering all possibilities, just the most likely input mistakes

    if islogical(ElementsToBeRefined)

        if numel(ElementsToBeRefined) ~= MUA.Nele

            fprintf("Number of elements in the logical list ElementsToBeRefined (%i) not equal to number of elements in mesh (%i) \n ",...
                numel(ElementsToBeRefined),MUA.Nele)
            error("Ua:IncorectUserInput","Incorrect DefineDesiredEleSize user input")
        end

    end

    if islogical(ElementsToBeCoarsened)

        if numel(ElementsToBeCoarsened) ~= MUA.Nele

            fprintf("Number of elements in the logical list ElementsToBeCoarsened (%i) not equal to number of elements in mesh (%i) \n ",...
                numel(ElementsToBeCoarsened),MUA.Nele)
            error("Ua:IncorectUserInput","Incorrect DefineDesiredEleSize user input")
        end

    end


elseif contains(CtrlVar.MeshRefinementMethod,'global','IgnoreCase',true)

    if ~isequal(InputElementsToBeRefined,ElementsToBeRefined)
        if nCalls==1
            fprintf(['Note: The mesh refinement methods is %s. \n' ...
                'When using this global mesh-refinement method in combination with ``DefineDesiredEleSize'' \n',...
                'you only need to specify as an output the variable ``EleSizeDesired''.\n',...
                'You can specify which elements are to be refined and/or coarsended but this will be ignored. \n'],...
                CtrlVar.MeshRefinementMethod)
        end
    end


else
    fprintf('Incorrect value for CtrlVar.MeshRefinementMethod (%s) \n',CtrlVar.MeshRefinementMethod)
    error('RemeshingBasedOnExplicitErrorEstimate:CaseNotFound','Case not found.')
end



end

